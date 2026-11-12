"""
Experimento controlado con una comunidad sintética donde SABEMOS la respuesta:
  X        lecturas de una especie inventada X (30 000 pb)
  X_rep    otras lecturas de la misma X (réplica técnica: mide el ruido de muestreo)
  X_cepa   lecturas de una cepa de X con una fracción d de SNPs
  Y        lecturas de otra especie, sin relación con X
Para varias divergencias d y varios k, calcula Bray-Curtis con Simka
(kmer_distances.py, idéntico a Simka) y con Nomka. Solo librería estándar.

Uso (desde la raíz del repositorio):
    python3 experimentos/sintetico.py
Salida:
    experimentos/tablas/sintetico.csv
"""

import os
import random
import sys

sys.path.insert(0, "src")
import kmer_distances as base
import Nomka

DIVERGENCIAS = [0.0, 0.005, 0.01, 0.02, 0.05]
VALORES_K = [15, 21, 31]
LARGO_GENOMA = 30000
COBERTURA = 15
LARGO_LECTURA = 101
T = 2


def genoma(largo):
    return "".join(random.choice("ACGT") for _ in range(largo))


def mutar(secuencia, tasa):
    """Cambia cada base por otra con probabilidad 'tasa' (SNPs)."""
    nueva = []
    for base_original in secuencia:
        if random.random() < tasa:
            nueva.append(random.choice([b for b in "ACGT" if b != base_original]))
        else:
            nueva.append(base_original)
    return "".join(nueva)


def lecturas(secuencia):
    """Lecturas al azar de cualquiera de las dos hebras, con la cobertura indicada."""
    resultado = []
    for _ in range(COBERTURA * len(secuencia) // LARGO_LECTURA):
        p = random.randint(0, len(secuencia) - LARGO_LECTURA)
        r = secuencia[p:p + LARGO_LECTURA]
        if random.random() < 0.5:
            r = base.reverse_complement(r)
        resultado.append(r)
    return resultado


def tablas_solidas(muestras, k):
    tablas = []
    for lista in muestras:
        conteos = {}
        for r in lista:
            base.add_read(r, k, conteos)
        base.remove_rare_kmers(conteos, T)
        tablas.append(conteos)
    return tablas


def main():
    random.seed(2026)
    os.makedirs("experimentos/tablas", exist_ok=True)
    x = genoma(LARGO_GENOMA)
    y = genoma(LARGO_GENOMA)
    lect_x, lect_rep, lect_y = lecturas(x), lecturas(x), lecturas(y)
    filas = ["divergencia,k,par,bc_simka,bc_nomka,jaccard_simka,jaccard_nomka"]
    for d in DIVERGENCIAS:
        lect_cepa = lecturas(mutar(x, d))
        muestras = [lect_x, lect_rep, lect_cepa, lect_y]
        for k in VALORES_K:
            tablas = tablas_solidas(muestras, k)
            bc_s, jac_s = base.distance_matrices(tablas)
            familias = Nomka.agrupar_variantes(tablas, Nomka.nombres_de_familia(tablas))
            bc_n, jac_n = base.distance_matrices(familias)
            for nombre, j in [("X-X_rep", 1), ("X-X_cepa", 2), ("X-Y", 3)]:
                filas.append("%g,%d,%s,%.6f,%.6f,%.6f,%.6f" % (d, k, nombre, bc_s[0][j], bc_n[0][j],
                                                              jac_s[0][j], jac_n[0][j]))
            print("d=%.3f k=%d  X-cepa: Simka %.3f Nomka %.3f | X-Y: Simka %.3f Nomka %.3f"
                  % (d, k, bc_s[0][2], bc_n[0][2], bc_s[0][3], bc_n[0][3]), flush=True)
    with open("experimentos/tablas/sintetico.csv", "w") as f:
        f.write("\n".join(filas) + "\n")


if __name__ == "__main__":
    main()
