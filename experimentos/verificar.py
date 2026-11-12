"""
Verificaciones experimentales, teóricas y computacionales de los resultados HMP.
Solo usa la librería estándar.

  1. Reproducibilidad: ¿las 3 repeticiones dan matrices idénticas (byte a byte)?
     Cada repetición es un proceso nuevo de Python, con una semilla de hash
     distinta (Python cambia la función hash de los strings en cada ejecución).
  2. Variabilidad de los recursos entre repeticiones (coeficiente de variación).
  3. Postcondiciones de la formalización: matriz simétrica, diagonal 0, valores en [0, 1].
  4. Propiedad teórica: Bray-Curtis de Nomka <= Bray-Curtis de Simka en toda pareja
     (min(a1 + a2, b1 + b2) >= min(a1, b1) + min(a2, b2) y los totales no cambian).
  5. Estabilidad frente a N: la distancia de una pareja no debería depender de qué
     otras muestras se incluyen. Se cumple en Simka; en Nomka puede cambiar (limitación).
  6. Implementación independiente: src/Simka.py (fases como el C++ original) debe dar
     lo mismo que kmer_distances.py sobre datos HMP reales (opción --independiente).

Uso (desde la raíz del repositorio):
    python3 experimentos/verificar.py [--independiente]
Salida:
    experimentos/tablas/verificacion.txt
"""

import hashlib
import os
import statistics
import sys

CARPETA = "resultados/hmp"
SALIDA = []


def anotar(texto=""):
    print(texto, flush=True)
    SALIDA.append(texto)


def leer_matriz(ruta):
    with open(ruta) as f:
        lineas = f.read().strip().split("\n")
    nombres = lineas[0].split(";")[1:]
    return nombres, [[float(x) for x in l.split(";")[1:]] for l in lineas[1:]]


def huella(ruta):
    with open(ruta, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def ruta(programa, k, n, rep, distancia):
    return "%s/%s_k%d_n%d_r%d_%s.csv" % (CARPETA, programa, k, n, rep, distancia)


def prueba_reproducibilidad():
    anotar("1. REPRODUCIBILIDAD DE LAS MATRICES (3 repeticiones, procesos y semillas de hash distintos)")
    for programa in ["Simka-Python", "Nomka"]:
        for k in [21, 31]:
            for n in [2, 4, 6, 9]:
                for distancia in ["braycurtis", "jaccard"]:
                    rutas = [ruta(programa, k, n, r, distancia) for r in [1, 2, 3]]
                    existentes = [r for r in rutas if os.path.exists(r)]
                    if len(existentes) < 2:
                        continue
                    distintas = len(set(huella(r) for r in existentes))
                    anotar("   %-12s k=%d n=%d %-10s repeticiones=%d  %s"
                           % (programa, k, n, distancia, len(existentes),
                              "IDÉNTICAS" if distintas == 1 else "¡DIFERENTES! (%d versiones)" % distintas))


def prueba_recursos():
    anotar("\n2. VARIABILIDAD DE TIEMPO Y RAM ENTRE REPETICIONES (CV = desviación / media)")
    with open("experimentos/tablas/mediciones.csv") as f:
        filas = [l.strip().split(",") for l in f.readlines()[1:]]
    for programa in ["Simka-Python", "Nomka"]:
        for n in ["2", "4", "6", "9"]:
            sel = [f for f in filas if f[0] == programa and f[1] == "21" and f[2] == n and f[10] == "0"]
            if len(sel) < 2:
                continue
            t = [float(f[4]) for f in sel]
            r = [float(f[7]) for f in sel]
            c = [float(f[6]) for f in sel]
            anotar("   %-12s n=%s  tiempo %.1f s (CV %.1f %%)  RAM %.0f MB (CV %.1f %%)  CPU media %.0f %%  [%d rep.]"
                   % (programa, n, statistics.mean(t), 100 * statistics.stdev(t) / statistics.mean(t),
                      statistics.mean(r), 100 * statistics.stdev(r) / statistics.mean(r),
                      statistics.mean(c), len(sel)))


def prueba_postcondiciones():
    anotar("\n3. POSTCONDICIONES (simétrica, diagonal 0, valores en [0, 1])")
    errores = 0
    revisadas = 0
    for nombre in sorted(os.listdir(CARPETA)):
        if not nombre.endswith(".csv"):
            continue
        nombres, d = leer_matriz(os.path.join(CARPETA, nombre))
        revisadas += 1
        for i in range(len(nombres)):
            if d[i][i] != 0:
                errores += 1
                anotar("   " + nombre + ": diagonal distinta de 0")
            for j in range(len(nombres)):
                if d[i][j] != d[j][i] or not (0 <= d[i][j] <= 1):
                    errores += 1
                    anotar("   %s: problema en (%s, %s) = %s" % (nombre, nombres[i], nombres[j], d[i][j]))
    anotar("   %d matrices revisadas, %d problemas" % (revisadas, errores))


def prueba_nomka_menor_o_igual():
    anotar("\n4. PROPIEDAD TEÓRICA: distancia Nomka <= distancia Simka en cada pareja")
    for k in [21, 31]:
        for n in [2, 4, 6, 9]:
            for distancia in ["braycurtis", "jaccard"]:
                rs, rn = ruta("Simka-Python", k, n, 1, distancia), ruta("Nomka", k, n, 1, distancia)
                if not (os.path.exists(rs) and os.path.exists(rn)):
                    continue
                nombres, ds = leer_matriz(rs)
                _, dn = leer_matriz(rn)
                pares = [(ds[i][j], dn[i][j]) for i in range(len(nombres)) for j in range(i + 1, len(nombres))]
                violaciones = sum(1 for s, m in pares if m > s + 1e-6)
                anotar("   k=%d n=%d %-10s parejas=%d  Nomka > Simka en %d  (Δ medio %.4f, Δ máx %.4f)"
                       % (k, n, distancia, len(pares), violaciones,
                          statistics.mean(s - m for s, m in pares), max(s - m for s, m in pares)))


def prueba_estabilidad_n():
    anotar("\n5. ESTABILIDAD FRENTE A N: ¿la distancia de una pareja cambia al agregar otras muestras?")
    for programa in ["Simka-Python", "Nomka"]:
        for distancia in ["braycurtis", "jaccard"]:
            valores = {}
            for n in [2, 4, 6, 9]:
                r = ruta(programa, 21, n, 1, distancia)
                if not os.path.exists(r):
                    continue
                nombres, d = leer_matriz(r)
                for i in range(len(nombres)):
                    for j in range(i + 1, len(nombres)):
                        valores.setdefault((nombres[i], nombres[j]), []).append(d[i][j])
            cambios = [max(v) - min(v) for v in valores.values() if len(v) > 1]
            if cambios:
                anotar("   %-12s %-10s parejas en varios N: %d  cambio máximo %.6f  cambio medio %.6f"
                       % (programa, distancia, len(cambios), max(cambios), statistics.mean(cambios)))


def prueba_independiente():
    anotar("\n6. IMPLEMENTACIÓN INDEPENDIENTE: src/Simka.py contra kmer_distances.py (HMP, n = 2, k = 21)")
    sys.path.insert(0, "src")
    import Simka
    muestras = {}
    with open("experimentos/listas/hmp_n2.txt") as f:
        for linea in f:
            nombre, archivo = [x.strip() for x in linea.split(":")]
            lecturas = []
            with open(os.path.join("experimentos/listas", archivo)) as fq:
                for p, l in enumerate(fq):
                    if p % 4 == 1:
                        lecturas.append(l.strip().upper())
            muestras[nombre] = lecturas
    resultado = Simka.simka(muestras, k=21, abundancia_min=2)
    for distancia, clave in [("braycurtis", "braycurtis"), ("jaccard", "jaccard")]:
        _, d = leer_matriz(ruta("Simka-Python", 21, 2, 1, distancia))
        anotar("   %-10s Simka.py = %.6f   kmer_distances.py = %.6f   %s"
               % (distancia, resultado[clave][0][1], d[0][1],
                  "IGUALES" if "%.6f" % resultado[clave][0][1] == "%.6f" % d[0][1] else "¡DISTINTOS!"))


def main():
    prueba_reproducibilidad()
    prueba_recursos()
    prueba_postcondiciones()
    prueba_nomka_menor_o_igual()
    prueba_estabilidad_n()
    if "--independiente" in sys.argv:
        prueba_independiente()
    os.makedirs("experimentos/tablas", exist_ok=True)
    with open("experimentos/tablas/verificacion.txt", "w") as f:
        f.write("\n".join(SALIDA) + "\n")


if __name__ == "__main__":
    main()
