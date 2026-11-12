"""
Compara las matrices de Simka (kmer_distances.py) y de Nomka sobre las muestras HMP.
Solo usa la librería estándar.

Métricas de calidad de cada matriz (el sitio corporal es el agrupamiento esperado):
  silueta           promedio de s(i) = (b - a) / max(a, b); a = distancia media a las
                    muestras de su mismo sitio, b = menor distancia media a otro sitio.
                    Cerca de 1 = los sitios quedan bien separados.
  vecino_ok         % de muestras cuya muestra más cercana es del mismo sitio
  dentro, entre     distancia media entre muestras del mismo sitio / de sitios distintos

Comparación Simka contra Nomka (por pareja de muestras):
  delta = distancia Simka - distancia Nomka (parte de la distancia que se debe a SNPs)
  spearman          correlación de rangos entre las dos matrices (1 = mismo orden)

Uso (desde la raíz del repositorio):
    python3 experimentos/comparar_matrices.py
Salida:
    experimentos/tablas/calidad.csv, experimentos/tablas/delta_pares.csv
"""

import os

CARPETA = "resultados/hmp"


def leer_matriz(ruta):
    """Lee una matriz con el formato de Simka (separador ';')."""
    with open(ruta) as f:
        lineas = f.read().strip().split("\n")
    nombres = lineas[0].split(";")[1:]
    matriz = []
    for linea in lineas[1:]:
        matriz.append([float(x) for x in linea.split(";")[1:]])
    return nombres, matriz


def sitio(nombre):
    """STOOL2 -> STOOL (se quitan los dígitos del final)."""
    return nombre.rstrip("0123456789")


def silueta(nombres, d):
    n = len(nombres)
    sitios = sorted(set(sitio(x) for x in nombres))
    total = 0.0
    for i in range(n):
        promedio = {}
        for s in sitios:
            otros = [d[i][j] for j in range(n) if j != i and sitio(nombres[j]) == s]
            if len(otros) > 0:
                promedio[s] = sum(otros) / len(otros)
        propio = sitio(nombres[i])
        if propio not in promedio:          # única muestra de su sitio: s(i) = 0
            continue
        a = promedio[propio]
        b = min(promedio[s] for s in promedio if s != propio)
        if max(a, b) > 0:
            total += (b - a) / max(a, b)
    return total / n


def vecino_correcto(nombres, d):
    n = len(nombres)
    aciertos = 0
    for i in range(n):
        mejor = None
        for j in range(n):
            if j != i and (mejor is None or d[i][j] < d[i][mejor]):
                mejor = j
        if sitio(nombres[mejor]) == sitio(nombres[i]):
            aciertos += 1
    return 100.0 * aciertos / n


def dentro_entre(nombres, d):
    dentro, entre = [], []
    for i in range(len(nombres)):
        for j in range(i + 1, len(nombres)):
            if sitio(nombres[i]) == sitio(nombres[j]):
                dentro.append(d[i][j])
            else:
                entre.append(d[i][j])
    media_dentro = sum(dentro) / len(dentro) if dentro else float("nan")
    media_entre = sum(entre) / len(entre) if entre else float("nan")
    return media_dentro, media_entre


def rangos(valores):
    """Rango de cada valor (empates reciben el rango promedio)."""
    orden = sorted(range(len(valores)), key=lambda i: valores[i])
    resultado = [0.0] * len(valores)
    i = 0
    while i < len(orden):
        j = i
        while j + 1 < len(orden) and valores[orden[j + 1]] == valores[orden[i]]:
            j += 1
        for p in range(i, j + 1):
            resultado[orden[p]] = (i + j) / 2 + 1
        i = j + 1
    return resultado


def spearman(x, y):
    rx, ry = rangos(x), rangos(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den > 0 else float("nan")


def pares(nombres, d):
    lista = []
    for i in range(len(nombres)):
        for j in range(i + 1, len(nombres)):
            lista.append((nombres[i], nombres[j], d[i][j]))
    return lista


def main():
    filas_calidad = ["k,n,distancia,programa,silueta,vecino_ok_pct,dentro,entre,spearman_vs_simka,delta_medio,delta_max"]
    filas_delta = ["k,distancia,muestra_a,muestra_b,mismo_sitio,simka,nomka,delta"]
    for k in [21, 31]:
        for n in [2, 4, 6, 9]:
            for distancia in ["braycurtis", "jaccard"]:
                ruta_s = CARPETA + "/Simka-Python_k%d_n%d_r1_%s.csv" % (k, n, distancia)
                ruta_n = CARPETA + "/Nomka_k%d_n%d_r1_%s.csv" % (k, n, distancia)
                if not (os.path.exists(ruta_s) and os.path.exists(ruta_n)):
                    continue
                nombres, ds = leer_matriz(ruta_s)
                _, dn = leer_matriz(ruta_n)
                ps = [p[2] for p in pares(nombres, ds)]
                pn = [p[2] for p in pares(nombres, dn)]
                deltas = [a - b for a, b in zip(ps, pn)]
                rho = spearman(ps, pn) if len(ps) > 2 else float("nan")
                for programa, d in [("Simka", ds), ("Nomka", dn)]:
                    sil = silueta(nombres, d)
                    vec = vecino_correcto(nombres, d)
                    den, ent = dentro_entre(nombres, d)
                    extra = ("1,0,0" if programa == "Simka"
                             else "%.4f,%.4f,%.4f" % (rho, sum(deltas) / len(deltas), max(deltas)))
                    filas_calidad.append("%d,%d,%s,%s,%.4f,%.1f,%.4f,%.4f,%s"
                                         % (k, n, distancia, programa, sil, vec, den, ent, extra))
                if n == 9:
                    for (a, b, vs), vn in zip(pares(nombres, ds), pn):
                        filas_delta.append("%d,%s,%s,%s,%s,%.6f,%.6f,%.6f"
                                           % (k, distancia, a, b, sitio(a) == sitio(b), vs, vn, vs - vn))
    os.makedirs("experimentos/tablas", exist_ok=True)
    with open("experimentos/tablas/calidad.csv", "w") as f:
        f.write("\n".join(filas_calidad) + "\n")
    with open("experimentos/tablas/delta_pares.csv", "w") as f:
        f.write("\n".join(filas_delta) + "\n")
    for fila in filas_calidad:
        print(fila.replace(",", "\t"))


if __name__ == "__main__":
    main()
