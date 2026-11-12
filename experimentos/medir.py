"""
Mide tiempo, CPU, RAM y disco de kmer_distances.py (mismos resultados que Simka)
y de Nomka.py sobre las muestras HMP. Solo usa la librería estándar.

Para cada corrida guarda:
  tiempo_s   tiempo de pared (reloj)
  cpu_s      tiempo de CPU (usuario + sistema) que reporta el sistema operativo
  cpu_pct    cpu_s / tiempo_s * 100 (100 % = un núcleo ocupado todo el tiempo)
  ram_mb     memoria RAM máxima del proceso (pico)
  entrada_mb tamaño de los archivos de lecturas leídos
  salida_mb  tamaño de las matrices escritas (no hay archivos temporales)

Uso (desde la raíz del repositorio):
    python3 experimentos/medir.py <corridas en paralelo>
Salida:
    experimentos/tablas/mediciones.csv
    resultados/hmp/<programa>_k<k>_n<n>_r<rep>_{braycurtis,jaccard}.csv
    resultados/hmp/logs/...
"""

import os
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

# Orden en que se agregan las muestras: se alternan los sitios corporales
ORDEN = ["STOOL1", "TONGUE1", "VAG1", "STOOL2", "TONGUE2", "VAG2", "STOOL3", "TONGUE3", "VAG3"]
TAMANOS = [2, 4, 6, 9]
PROGRAMAS = {"Simka-Python": "src/kmer_distances.py", "Nomka": "src/Nomka.py"}
T = 2
SALIDA = "resultados/hmp"
TABLA = "experimentos/tablas/mediciones.csv"

# Como mucho 2 corridas grandes de Nomka a la vez, para no llenar la RAM
SEMAFORO_GRANDES = threading.Semaphore(2)
CANDADO_TABLA = threading.Lock()


def escribir_listas():
    """Crea experimentos/listas/hmp_n<n>.txt con las primeras n muestras de ORDEN."""
    os.makedirs("experimentos/listas", exist_ok=True)
    for n in TAMANOS:
        with open("experimentos/listas/hmp_n" + str(n) + ".txt", "w") as f:
            for nombre in ORDEN[:n]:
                f.write(nombre + ": ../../datos_hmp/" + nombre + ".fastq\n")


def megabytes(rutas):
    total = 0
    for ruta in rutas:
        if os.path.exists(ruta):
            total += os.path.getsize(ruta)
    return total / 1048576


def correr(programa, k, n, rep):
    lista = "experimentos/listas/hmp_n" + str(n) + ".txt"
    prefijo = SALIDA + "/" + programa + "_k" + str(k) + "_n" + str(n) + "_r" + str(rep)
    log = SALIDA + "/logs/" + programa + "_k" + str(k) + "_n" + str(n) + "_r" + str(rep) + ".log"
    grande = programa == "Nomka" and n >= 6
    if grande:
        SEMAFORO_GRANDES.acquire()
    try:
        inicio = time.time()
        with open(log, "w") as salida_log:
            proceso = subprocess.Popen(["python3", PROGRAMAS[programa], lista, str(k), str(T), prefijo],
                                       stdout=salida_log, stderr=subprocess.STDOUT)
            _, estado, uso = os.wait4(proceso.pid, 0)
        fin = time.time()
    finally:
        if grande:
            SEMAFORO_GRANDES.release()

    tiempo = fin - inicio
    cpu = uso.ru_utime + uso.ru_stime
    ram_mb = uso.ru_maxrss / 1048576          # en macOS ru_maxrss viene en bytes
    entrada = megabytes(["datos_hmp/" + m + ".fastq" for m in ORDEN[:n]])
    salida = megabytes([prefijo + "_braycurtis.csv", prefijo + "_jaccard.csv"])
    fila = [programa, k, n, rep, "%.2f" % tiempo, "%.2f" % cpu, "%.1f" % (100 * cpu / tiempo),
            "%.1f" % ram_mb, "%.1f" % entrada, "%.4f" % salida, "%d" % os.waitstatus_to_exitcode(estado),
            time.strftime("%H:%M:%S", time.localtime(inicio)), time.strftime("%H:%M:%S", time.localtime(fin))]
    with CANDADO_TABLA:
        with open(TABLA, "a") as f:
            f.write(",".join(str(x) for x in fila) + "\n")
    print(" ".join(str(x) for x in fila), flush=True)


def main():
    global TABLA
    en_paralelo = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    control = "control" in sys.argv
    if control:
        # Corridas de control: cada programa SOLO (sin otras corridas en paralelo),
        # para comprobar que el paralelismo no sesgó el tiempo ni la RAM
        TABLA = "experimentos/tablas/mediciones_control.csv"
        en_paralelo = 1
    escribir_listas()
    os.makedirs(SALIDA + "/logs", exist_ok=True)
    os.makedirs("experimentos/tablas", exist_ok=True)
    if not os.path.exists(TABLA):
        with open(TABLA, "w") as f:
            f.write("programa,k,n,rep,tiempo_s,cpu_s,cpu_pct,ram_mb,entrada_mb,salida_mb,codigo_salida,inicio,fin\n")

    # Plan: k = 21 con n = 2, 4, 6, 9 y 3 repeticiones; k = 31 solo con n = 9 y 1 repetición
    corridas = []
    for rep in [1, 2, 3]:
        for n in TAMANOS:
            for programa in PROGRAMAS:
                corridas.append((programa, 21, n, rep))
    for programa in PROGRAMAS:
        corridas.append((programa, 31, 9, 1))
    # Primero la repetición 1 completa (así hay una serie completa cuanto antes),
    # luego las repeticiones 2 y 3, y al final k = 31. Dentro de cada grupo,
    # las corridas más largas primero.
    def prioridad(c):
        programa, k, n, rep = c
        grupo = 99 if k == 31 else rep
        return (grupo, -n, programa != "Nomka")
    corridas.sort(key=prioridad)
    if control:
        corridas = [("Simka-Python", 21, 9, "C"), ("Nomka", 21, 9, "C")]

    print("Corridas:", len(corridas), " en paralelo:", en_paralelo, flush=True)
    with ThreadPoolExecutor(max_workers=en_paralelo) as grupo:
        for c in corridas:
            grupo.submit(correr, *c)


if __name__ == "__main__":
    main()
