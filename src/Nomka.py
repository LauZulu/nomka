"""
Nomka: compara N muestras metagenómicas como Simka, pero tolerando SNPs.

Simka (y kmer_distances.py) solo considera "compartidas" dos palabras de ADN
(k-mers) si son idénticas letra por letra. Nomka agrega UN paso antes de
comparar: junta en una "familia" los k-mers que difieren en una sola letra,
porque lo más probable es que sean la misma región de la misma especie en
cepas distintas (un SNP). Todo lo demás es igual: k-mers canónicos, filtro t,
Bray-Curtis y Jaccard.

Uso (desde la raíz del repositorio):
    python3 src/Nomka.py <archivo de muestras> <k> <t> <prefijo de salida>

Ejemplo:
    python3 src/Nomka.py datos_simka/simka_input.txt 21 2 resultados/nomka_ejemplo

El archivo de muestras tiene el mismo formato que Simka:
una muestra por línea, "nombre: archivo1 , archivo2 ; archivo3".
"""

import sys

# Se reutiliza el código de Laura sin cambios: leer, contar, filtrar, comparar y escribir.
from kmer_distances import (read_samples_file, count_file, remove_rare_kmers,
                            canonical, distance_matrices, write_matrix)


# ---------------------------------------------------------------------------
# Lo nuevo de Nomka: agrupar los k-mers que difieren en un SNP
# ---------------------------------------------------------------------------

def variantes(kmer):
    """
    Genera los "primos" de un k-mer: todas las palabras del mismo largo que
    difieren en UNA sola letra. Con k = 21 son 21 posiciones x 3 letras = 63.
    Ejemplo (k = 3): ACG -> CCG, GCG, TCG, AAG, AGG, ATG, ACA, ACC, ACT
    :param kmer: k-mer canónico
    :return: lista con la forma canónica de cada primo
    """
    lista = []
    for posicion in range(len(kmer)):
        for letra in "ACGT":
            if letra != kmer[posicion]:
                primo = kmer[:posicion] + letra + kmer[posicion + 1:]
                # El primo también se pasa a canónico, para reconocerlo aunque
                # se haya leído desde la otra hebra
                lista.append(canonical(primo))
    return lista


def nombres_de_familia(tablas):
    """
    Decide a qué familia pertenece cada k-mer sólido.
    El nombre de la familia es el menor (en orden alfabético) entre el k-mer y
    sus primos que existen en alguna muestra. Es una regla fija, para que
    todas las muestras usen la misma etiqueta para la misma familia.
    :param tablas: tablas[i] tiene los conteos de k-mers sólidos de la muestra i
    :return: diccionario k-mer -> nombre de su familia
    """
    # 1. Todos los k-mers sólidos de todas las muestras
    todos = {}
    for tabla in tablas:
        for kmer in tabla:
            todos[kmer] = True

    # 2. Para cada k-mer, buscar cuáles de sus primos existen y quedarse con el menor
    nombre = {}
    for kmer in todos:
        menor = kmer
        for primo in variantes(kmer):
            if primo in todos and primo < menor:
                menor = primo
        nombre[kmer] = menor
    return nombre


def agrupar_variantes(tablas, nombre):
    """
    Reescribe cada tabla usando el nombre de familia y suma los conteos.
    Ejemplo: {ACGTA: 6, TTGCA: 3} con nombre[ACGTA] = ACCTA -> {ACCTA: 6, TTGCA: 3}
    :param tablas: tablas[i] tiene los conteos de k-mers sólidos de la muestra i
    :param nombre: diccionario k-mer -> nombre de su familia
    :return: tablas nuevas, con familias en vez de k-mers
    """
    nuevas = []
    for tabla in tablas:
        nueva = {}
        for kmer in tabla:
            familia = nombre[kmer]
            nueva[familia] = nueva.get(familia, 0) + tabla[kmer]
        nuevas.append(nueva)
    return nuevas


def resumen_de_familias(nombre):
    """
    Cuenta cuántos k-mers quedaron en cada familia, para revisar el resultado.
    :param nombre: diccionario k-mer -> nombre de su familia
    :return: (número de familias, familias con más de un k-mer, tamaño de la familia más grande)
    """
    tamano = {}
    for kmer in nombre:
        familia = nombre[kmer]
        tamano[familia] = tamano.get(familia, 0) + 1
    con_variantes = 0
    mas_grande = 0
    for familia in tamano:
        if tamano[familia] > 1:
            con_variantes += 1
        if tamano[familia] > mas_grande:
            mas_grande = tamano[familia]
    return len(tamano), con_variantes, mas_grande


# ---------------------------------------------------------------------------
# Programa principal
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 5:
        print("Uso: python3 src/Nomka.py <archivo de muestras> <k> <t> <prefijo de salida>")
        sys.exit(1)
    archivo_muestras = sys.argv[1]
    k = int(sys.argv[2])
    t = int(sys.argv[3])
    prefijo = sys.argv[4]

    # 1. Leer la lista de muestras (igual que kmer_distances.py)
    nombres, archivos = read_samples_file(archivo_muestras)

    # 2. Contar los k-mers de cada muestra y quedarse con los sólidos (igual que kmer_distances.py)
    tablas = []
    for i in range(len(nombres)):
        conteos = {}
        for archivo in archivos[i]:
            count_file(archivo, k, conteos)
        remove_rare_kmers(conteos, t)
        tablas.append(conteos)

    # 3. NUEVO: agrupar en familias los k-mers que difieren en un SNP
    nombre = nombres_de_familia(tablas)
    familias = agrupar_variantes(tablas, nombre)

    num_familias, con_variantes, mas_grande = resumen_de_familias(nombre)
    print("k-mers sólidos distintos (todas las muestras): " + str(len(nombre)))
    print("Familias: " + str(num_familias) + " (con variantes: " + str(con_variantes)
          + "; la más grande tiene " + str(mas_grande) + " k-mers)")
    print("Muestra\tKmersSolidos\tFamilias\tTotal")
    for i in range(len(nombres)):
        print(nombres[i] + "\t" + str(len(tablas[i])) + "\t" + str(len(familias[i]))
              + "\t" + str(sum(familias[i].values())))

    # 4. Bray-Curtis y Jaccard sobre las familias (igual que kmer_distances.py)
    bray_curtis, jaccard = distance_matrices(familias)
    write_matrix(prefijo + "_braycurtis.csv", nombres, bray_curtis)
    write_matrix(prefijo + "_jaccard.csv", nombres, jaccard)
    print("Matrices escritas con el prefijo " + prefijo)


if __name__ == "__main__":
    main()
