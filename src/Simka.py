"""
Simka en Python (versión de aprendizaje).

Sigue la lógica del repositorio oficial github.com/GATB/simka, sin la parte
de ingeniería (particiones en disco, paralelismo, compresión).

Fases (y dónde están en el repositorio original):
    0. Leer las muestras y filtrar lecturas  -> SimkaCommons.hpp (SimkaSequenceFilter)
    1. Contar k-mers de cada muestra          -> SimkaCount.cpp / MiniKC.hpp
       y quedarse con los k-mers "sólidos"
    2. Unir los conteos en un vector por      -> SimkaMerge.cpp / SimkaAlgorithm.hpp
       k-mer y acumular estadísticas por        (updateDistanceDefault)
       pareja de muestras
    3. Calcular las distancias al final       -> SimkaDistance.cpp
"""

import math

# Valores por defecto, iguales a los de Simka (src/core/Simka.cpp)
K_POR_DEFECTO = 21
ABUNDANCIA_MIN_POR_DEFECTO = 2
ABUNDANCIA_MAX_POR_DEFECTO = 999999999


# ===========================================================================
# FASE 0: leer muestras y filtrar lecturas
# ===========================================================================

def leer_fasta(ruta: str) -> list:
    """
    Lee un archivo FASTA y retorna la lista de secuencias (lecturas).
    Las líneas que empiezan con '>' son encabezados; las demás son ADN.
    Una secuencia puede ocupar varias líneas.
    """
    lecturas = []
    actual = ""
    archivo = open(ruta, "r")
    for linea in archivo:
        linea = linea.strip()
        if linea.startswith(">"):
            if actual != "":
                lecturas.append(actual)
            actual = ""
        else:
            actual += linea.upper()
    if actual != "":
        lecturas.append(actual)
    archivo.close()
    return lecturas


def leer_archivo_entrada(ruta: str) -> dict:
    """
    Lee el archivo de entrada con el formato de Simka:
        A: A.fasta
        D: D_paired_1.fasta ; D_paired_2.fasta
        E: A.fasta , A.fasta ; B.fasta , B.fasta
    ',' separa partes de una muestra y ';' separa archivos pareados.
    En ambos casos Simka une todas las lecturas en una sola muestra.

    Retorna {id_muestra: [lecturas]}.
    """
    carpeta = ""
    if "/" in ruta:
        carpeta = ruta[:ruta.rfind("/") + 1]

    muestras = {}
    archivo = open(ruta, "r")
    for linea in archivo:
        linea = linea.strip()
        if linea == "":
            continue
        partes = linea.split(":")
        id_muestra = partes[0].strip()
        nombres_archivos = partes[1].replace(";", ",").split(",")

        lecturas = []
        for nombre in nombres_archivos:
            lecturas += leer_fasta(carpeta + nombre.strip())
        muestras[id_muestra] = lecturas
    archivo.close()
    return muestras


def indice_shannon(secuencia: str) -> float:
    """
    Mide qué tan "variada" es una secuencia (0 = una sola letra, 2 = A, C, G
    y T en igual proporción). Sirve para descartar lecturas de baja
    complejidad como "AAAAAAAAAA".
    """
    frecuencias = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
    for letra in secuencia:
        if letra in frecuencias:
            frecuencias[letra] += 1
        else:
            frecuencias["N"] += 1

    indice = 0.0
    for letra in frecuencias:
        p = frecuencias[letra] / len(secuencia)
        if p != 0:
            indice += p * math.log2(p)
    return abs(indice)


def filtrar_lecturas(lecturas: list, tamano_min: int, shannon_min: float) -> list:
    """
    Igual que SimkaSequenceFilter: descarta lecturas más cortas que
    `tamano_min` o con índice de Shannon menor que `shannon_min`.
    Con valor 0 el filtro está apagado (es el default en Simka).
    """
    rta = []
    for lectura in lecturas:
        if tamano_min != 0 and len(lectura) < tamano_min:
            continue
        if shannon_min != 0 and indice_shannon(lectura) < shannon_min:
            continue
        rta.append(lectura)
    return rta


# ===========================================================================
# FASE 1: contar k-mers de cada muestra (SimkaCount)
# ===========================================================================

def reverso_complementario(secuencia: str) -> str:
    """
    Devuelve la hebra complementaria leída al revés.
    Ejemplo: "AACG" -> complemento "TTGC" -> al revés "CGTT".
    """
    complemento = {"A": "T", "T": "A", "C": "G", "G": "C"}
    rta = ""
    for i in range(len(secuencia) - 1, -1, -1):
        rta += complemento[secuencia[i]]
    return rta


def kmer_canonico(kmer: str) -> str:
    """
    Un k-mer y su reverso complementario son el mismo trozo de ADN visto
    desde cada hebra. Simka los cuenta como uno solo (ModelCanonical de
    GATB), quedándose con el "menor" de los dos.
    """
    rc = reverso_complementario(kmer)
    if rc < kmer:
        return rc
    return kmer


def es_kmer_valido(kmer: str) -> bool:
    """Solo aceptamos k-mers hechos de A, C, G, T (simplificación nuestra)."""
    for letra in kmer:
        if letra not in "ACGT":
            return False
    return True


def contar_kmers(lecturas: list, k: int) -> dict:
    """
    Recorre cada lectura con una ventana de tamaño k y cuenta cuántas veces
    aparece cada k-mer canónico. Retorna {kmer: conteo}.
    """
    conteos = {}
    for lectura in lecturas:
        for i in range(len(lectura) - k + 1):
            kmer = lectura[i:i + k]
            if not es_kmer_valido(kmer):
                continue
            kmer = kmer_canonico(kmer)
            if kmer in conteos:
                conteos[kmer] += 1
            else:
                conteos[kmer] = 1
    return conteos


def filtrar_kmers_solidos(conteos: dict, abundancia_min: int, abundancia_max: int) -> dict:
    """
    Un k-mer es "sólido" en una muestra si su conteo está en
    [abundancia_min, abundancia_max] (MiniKC.hpp, línea 56).
    - Muy pocos -> probable error de secuenciación.
    - Demasiados -> probable contaminante o repetición.
    """
    rta = {}
    for kmer in conteos:
        if abundancia_min <= conteos[kmer] <= abundancia_max:
            rta[kmer] = conteos[kmer]
    return rta


def totales_por_muestra(solidos: dict) -> tuple:
    """
    Simka guarda dos números por muestra mientras cuenta:
      - total_kmers:     suma de las abundancias      (_nbSolidKmersPerBank)
      - kmers_distintos: cuántos k-mers diferentes hay (_nbSolidDistinctKmersPerBank)
    """
    total_kmers = 0
    for kmer in solidos:
        total_kmers += solidos[kmer]
    kmers_distintos = len(solidos)
    return total_kmers, kmers_distintos


# ===========================================================================
# FASE 2: unir conteos y acumular estadísticas por pareja (SimkaMerge)
# ===========================================================================

def construir_vectores_de_conteo(solidos_por_muestra: list) -> dict:
    """
    Une los conteos de todas las muestras: para cada k-mer, un vector con
    su abundancia en cada muestra (0 si no está). Es la "matriz de
    abundancias" del artículo: fila = k-mer, columna = muestra.
    """
    n = len(solidos_por_muestra)
    vectores = {}
    for j in range(n):
        for kmer in solidos_por_muestra[j]:
            if kmer not in vectores:
                vectores[kmer] = [0] * n
            vectores[kmer][j] = solidos_por_muestra[j][kmer]
    return vectores


def matriz_ceros(n: int) -> list:
    matriz = []
    for _ in range(n):
        matriz.append([0] * n)
    return matriz


def acumular_estadisticas(vectores: dict, n: int) -> tuple:
    """
    La idea central de Simka: se recorre la matriz UNA sola vez, k-mer por
    k-mer, y con cada vector se actualizan TODAS las parejas de muestras.

    Por cada k-mer solo importan las muestras donde aparece (conteo > 0);
    las demás parejas no cambian (updateDistanceDefault).

    Acumula, para cada pareja (i, j):
      - compartidos[i][j]: cuántos k-mers distintos comparten  (_matrixNbDistinctSharedKmers)
      - suma_min[i][j]:    suma de min(conteo_i, conteo_j)       (_brayCurtisNumerator)
    """
    compartidos = matriz_ceros(n)
    suma_min = matriz_ceros(n)

    for kmer in vectores:
        conteos = vectores[kmer]

        # Muestras donde aparece este k-mer ("_sharedBanks")
        presentes = []
        for i in range(n):
            if conteos[i] > 0:
                presentes.append(i)

        # Solo las parejas que comparten el k-mer
        for a in range(len(presentes)):
            for b in range(a + 1, len(presentes)):
                i = presentes[a]
                j = presentes[b]
                compartidos[i][j] += 1
                if conteos[i] < conteos[j]:
                    suma_min[i][j] += conteos[i]
                else:
                    suma_min[i][j] += conteos[j]

    return compartidos, suma_min


# ===========================================================================
# FASE 3: distancias (SimkaDistance)
# ===========================================================================

def distancia_bray_curtis(suma_min_ij: int, total_i: int, total_j: int) -> float:
    """
    BC = 1 - 2 * suma(min(a, b)) / (total_i + total_j)
    Usa abundancias. Si ambas muestras están vacías, Simka devuelve 1.
    """
    if total_i + total_j == 0:
        return 1.0
    return 1 - (2 * suma_min_ij) / (total_i + total_j)


def distancia_jaccard(compartidos_ij: int, distintos_i: int, distintos_j: int) -> float:
    """
    Simka la escribe con a, b, c (get_abc + jaccardCanberra):
      a = k-mers en ambas muestras
      b = k-mers solo en i  = distintos_i - a
      c = k-mers solo en j  = distintos_j - a
      J = (b + c) / (a + b + c)       <- equivale a 1 - |A∩B| / |A∪B|
    Solo mira presencia/ausencia.
    """
    a = compartidos_ij
    b = distintos_i - a
    c = distintos_j - a
    if a + b + c == 0:
        return 1.0
    return (b + c) / (a + b + c)


def simka(muestras: dict,
          k: int = K_POR_DEFECTO,
          abundancia_min: int = ABUNDANCIA_MIN_POR_DEFECTO,
          abundancia_max: int = ABUNDANCIA_MAX_POR_DEFECTO,
          tamano_min_lectura: int = 0,
          shannon_min_lectura: float = 0) -> dict:
    """
    Ejecuta Simka completo.
    muestras: {id_muestra: [lecturas]}
    Retorna {"nombres": [...], "braycurtis": matriz, "jaccard": matriz}
    """
    nombres = list(muestras.keys())
    n = len(nombres)

    # Fases 0 y 1: por cada muestra, filtrar lecturas, contar y filtrar k-mers
    solidos_por_muestra = []
    total_kmers = []
    kmers_distintos = []
    for nombre in nombres:
        lecturas = filtrar_lecturas(muestras[nombre], tamano_min_lectura, shannon_min_lectura)
        conteos = contar_kmers(lecturas, k)
        solidos = filtrar_kmers_solidos(conteos, abundancia_min, abundancia_max)
        total, distintos = totales_por_muestra(solidos)
        solidos_por_muestra.append(solidos)
        total_kmers.append(total)
        kmers_distintos.append(distintos)

    # Fase 2: una sola pasada por los k-mers
    vectores = construir_vectores_de_conteo(solidos_por_muestra)
    compartidos, suma_min = acumular_estadisticas(vectores, n)

    # Fase 3: distancias (la diagonal queda en 0)
    braycurtis = matriz_ceros(n)
    jaccard = matriz_ceros(n)
    for i in range(n):
        for j in range(i + 1, n):
            bc = distancia_bray_curtis(suma_min[i][j], total_kmers[i], total_kmers[j])
            jc = distancia_jaccard(compartidos[i][j], kmers_distintos[i], kmers_distintos[j])
            braycurtis[i][j] = bc
            braycurtis[j][i] = bc
            jaccard[i][j] = jc
            jaccard[j][i] = jc

    return {"nombres": nombres, "braycurtis": braycurtis, "jaccard": jaccard}


# ===========================================================================
# Salida: mismo formato que los .csv de Simka (separador ';', 6 decimales)
# ===========================================================================

def matriz_a_texto(nombres: list, matriz: list) -> str:
    texto = ""
    for nombre in nombres:
        texto += ";" + nombre
    texto += "\n"
    for i in range(len(nombres)):
        texto += nombres[i]
        for j in range(len(nombres)):
            texto += ";" + "%f" % matriz[i][j]
        texto += "\n"
    return texto


def comparar_con_referencia(texto: str, ruta_referencia: str) -> bool:
    """Compara nuestra matriz con la que produjo el Simka original."""
    archivo = open(ruta_referencia, "r")
    referencia = archivo.read()
    archivo.close()
    return texto == referencia


# ===========================================================================
# Ejemplo: los mismos datos y parámetros que tests/simple_test.py del repo
# ===========================================================================

if __name__ == "__main__":
    muestras = leer_archivo_entrada("datos_simka/simka_input.txt")

    for abundancia_min in [2, 0]:
        print("=== k = 21, abundancia_min =", abundancia_min, "===")
        resultado = simka(muestras, k=21, abundancia_min=abundancia_min)
        carpeta_ref = "datos_simka/referencia_k21_t" + str(abundancia_min) + "/"

        for distancia, archivo in [("braycurtis", "mat_abundance_braycurtis.csv"),
                                   ("jaccard", "mat_presenceAbsence_jaccard.csv")]:
            texto = matriz_a_texto(resultado["nombres"], resultado[distancia])
            igual = comparar_con_referencia(texto, carpeta_ref + archivo)
            print("\n" + archivo, "-> igual al Simka original:", igual)
            print(texto)
