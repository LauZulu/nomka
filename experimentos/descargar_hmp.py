"""
Descarga las primeras N lecturas (archivo R1) de las 9 muestras HMP del proyecto.
Python coordina la descarga y llama a fastq-dump (sra-tools, `brew install
sratoolkit`), que trae solo las primeras N lecturas de cada run desde NCBI, sin
bajar los archivos completos. (ENA no tiene el FASTQ de estos runs.)

Uso (desde la raíz del repositorio):
    python3 experimentos/descargar_hmp.py 100000
Salida:
    datos_hmp/<muestra>.fastq   (la carpeta no se versiona en git)
"""

import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

# Las 9 muestras HMP (proyecto SRA SRP002163), 3 por sitio corporal
MUESTRAS = [
    ("STOOL1", "SRR1031071"), ("STOOL2", "SRR1031141"), ("STOOL3", "SRR1031252"),
    ("TONGUE1", "SRR1031078"), ("TONGUE2", "SRR1031268"), ("TONGUE3", "SRR1031272"),
    ("VAG1", "SRR1031066"), ("VAG2", "SRR1564168"), ("VAG3", "SRR1565738"),
]
CARPETA = "datos_hmp"


def descargar(nombre, run, n_lecturas):
    destino = os.path.join(CARPETA, nombre + ".fastq")
    if os.path.exists(destino):
        return nombre + ": ya existe"
    temporal = os.path.join(CARPETA, "tmp_" + nombre)
    subprocess.run(["fastq-dump", "-X", str(n_lecturas), "--split-files", "-O", temporal, run],
                   check=True, capture_output=True)
    # Solo se usa el archivo R1 (como en el plan de la Entrega 3)
    os.rename(os.path.join(temporal, run + "_1.fastq"), destino)
    for archivo in os.listdir(temporal):
        os.remove(os.path.join(temporal, archivo))
    os.rmdir(temporal)
    lecturas = 0
    with open(destino) as f:
        for _ in f:
            lecturas += 1
    return nombre + " (" + run + "): " + str(lecturas // 4) + " lecturas"


def main():
    n_lecturas = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
    os.makedirs(CARPETA, exist_ok=True)
    with ThreadPoolExecutor(max_workers=9) as grupo:
        tareas = [grupo.submit(descargar, nombre, run, n_lecturas) for nombre, run in MUESTRAS]
        for tarea in tareas:
            print(tarea.result(), flush=True)


if __name__ == "__main__":
    main()
