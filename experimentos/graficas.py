"""
Gráficas y tabla resumen para el informe (etiquetas en inglés).
Usa matplotlib y scipy (solo para dibujar), instalados en .venv:
    python3 -m venv .venv && .venv/bin/pip install matplotlib scipy
Uso (desde la raíz del repositorio):
    .venv/bin/python experimentos/graficas.py
Lee experimentos/tablas/*.csv y resultados/hmp/*.csv; escribe experimentos/figuras/*.png
y experimentos/tablas/resumen.md
"""

import csv
import os
import statistics

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

FIG = "experimentos/figuras"
COLOR = {"Simka-Python": "#1f77b4", "Simka": "#1f77b4", "Nomka": "#d62728"}
ETIQUETA = {"Simka-Python": "Simka (Python, kmer_distances.py)", "Nomka": "Nomka"}


def leer_csv(ruta):
    with open(ruta) as f:
        return list(csv.DictReader(f))


def leer_matriz(ruta):
    with open(ruta) as f:
        lineas = f.read().strip().split("\n")
    nombres = lineas[0].split(";")[1:]
    return nombres, [[float(x) for x in l.split(";")[1:]] for l in lineas[1:]]


def media_de(valores):
    return statistics.mean(valores), (statistics.stdev(valores) if len(valores) > 1 else 0.0)


def figura_recursos(med):
    fig, ejes = plt.subplots(1, 3, figsize=(13, 3.8))
    for programa in ["Simka-Python", "Nomka"]:
        ns, t, te, r, re, c, ce = [], [], [], [], [], [], []
        for n in [2, 4, 6, 9]:
            filas = [f for f in med if f["programa"] == programa and f["k"] == "21" and f["n"] == str(n)]
            if not filas:
                continue
            ns.append(n)
            m, s = media_de([float(f["tiempo_s"]) for f in filas]); t.append(m); te.append(s)
            m, s = media_de([float(f["ram_mb"]) for f in filas]); r.append(m); re.append(s)
            m, s = media_de([float(f["cpu_s"]) for f in filas]); c.append(m); ce.append(s)
        for eje, y, e in [(ejes[0], t, te), (ejes[1], c, ce), (ejes[2], r, re)]:
            eje.errorbar(ns, y, yerr=e, marker="o", capsize=4, color=COLOR[programa], label=ETIQUETA[programa])
    for eje, titulo in zip(ejes, ["Wall-clock time (s)", "CPU time (s)", "Peak RAM (MB)"]):
        eje.set_xlabel("Number of samples (N)")
        eje.set_title(titulo)
        eje.set_xticks([2, 4, 6, 9])
        eje.grid(alpha=0.3)
    ejes[0].legend()
    fig.suptitle("Computational cost, k = 21, 100,000 reads per sample (mean ± SD of 3 runs)")
    fig.tight_layout()
    fig.savefig(FIG + "/fig_recursos.png", dpi=200)


def figura_matrices(k=21, distancia="braycurtis"):
    fig = plt.figure(figsize=(12, 7))
    for col, programa in enumerate(["Simka-Python", "Nomka"]):
        ruta = "resultados/hmp/%s_k%d_n9_r1_%s.csv" % (programa, k, distancia)
        if not os.path.exists(ruta):
            continue
        nombres, d = leer_matriz(ruta)
        z = linkage(squareform(d, checks=False), method="average")   # 'average' = UPGMA
        eje_arbol = fig.add_axes([0.08 + col * 0.47, 0.74, 0.34, 0.2])
        orden = dendrogram(z, labels=nombres, ax=eje_arbol, color_threshold=0, above_threshold_color="k",
                           no_labels=True)["leaves"]
        eje_arbol.axis("off")
        eje_arbol.set_title(ETIQUETA[programa] + " — UPGMA")
        eje = fig.add_axes([0.08 + col * 0.47, 0.17, 0.34, 0.56])
        m = [[d[i][j] for j in orden] for i in orden]
        im = eje.imshow(m, cmap="viridis_r", vmin=0, vmax=1)
        eje.set_xticks(range(len(orden)))
        eje.set_xticklabels([nombres[i] for i in orden], rotation=90, fontsize=8)
        eje.set_yticks(range(len(orden)))
        eje.set_yticklabels([nombres[i] for i in orden], fontsize=8)
        for i in range(len(orden)):
            for j in range(len(orden)):
                eje.text(j, i, "%.2f" % m[i][j], ha="center", va="center", fontsize=6,
                         color="white" if m[i][j] < 0.5 else "black")
    barra = fig.add_axes([0.93, 0.17, 0.015, 0.56])
    fig.colorbar(im, cax=barra, label="Bray-Curtis dissimilarity")
    fig.savefig(FIG + "/fig_matrices_k%d_%s.png" % (k, distancia), dpi=200)


def figura_delta(k=21):
    """Δ = Simka - Nomka por pareja, separado en mismo sitio / sitios distintos."""
    import random
    random.seed(1)
    filas = [f for f in leer_csv("experimentos/tablas/delta_pares.csv") if f["k"] == str(k)]
    fig, ejes = plt.subplots(1, 2, figsize=(10, 4.2))
    for eje, distancia, titulo in [(ejes[0], "braycurtis", "Bray-Curtis"), (ejes[1], "jaccard", "Jaccard")]:
        for x, (mismo, color, etiqueta) in enumerate([("True", "#2ca02c", "Same body site"),
                                                      ("False", "#7f7f7f", "Different body site")]):
            valores = [float(f["delta"]) for f in filas if f["distancia"] == distancia and f["mismo_sitio"] == mismo]
            eje.scatter([x + random.uniform(-0.12, 0.12) for _ in valores], valores, color=color, s=30, alpha=0.8)
            media = statistics.mean(valores)
            eje.plot([x - 0.25, x + 0.25], [media, media], color="k", lw=2)
            eje.text(x + 0.28, media, "mean %.4f\n(%d pairs)" % (media, len(valores)), va="center", fontsize=8)
        eje.axhline(0, color="k", lw=0.8, ls="--")
        eje.set_xticks([0, 1])
        eje.set_xticklabels(["Same body site", "Different body site"])
        eje.set_xlim(-0.5, 1.9)
        eje.set_ylabel("Δ = Simka − Nomka")
        eje.set_title("%s, k = %d, N = 9 (36 pairs)" % (titulo, k))
        eje.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(FIG + "/fig_delta_k%d.png" % k, dpi=200)


def figura_sintetico():
    filas = leer_csv("experimentos/tablas/sintetico.csv")
    fig, ejes = plt.subplots(1, 2, figsize=(11, 4.2))
    for k, estilo in [(15, ":"), (21, "-"), (31, "--")]:
        sel = [f for f in filas if f["k"] == str(k) and f["par"] == "X-X_cepa"]
        d = [100 * float(f["divergencia"]) for f in sel]
        ejes[0].plot(d, [float(f["bc_simka"]) for f in sel], estilo, marker="o", color=COLOR["Simka"],
                     label="Simka, k = %d" % k)
        ejes[0].plot(d, [float(f["bc_nomka"]) for f in sel], estilo, marker="s", color=COLOR["Nomka"],
                     label="Nomka, k = %d" % k)
    piso = [float(f["bc_simka"]) for f in filas if f["k"] == "21" and f["par"] == "X-X_rep"]
    ejes[0].axhline(statistics.mean(piso), color="gray", lw=1, label="Technical replicate (noise floor)")
    ejes[0].set_xlabel("SNP divergence between strains X and X' (%)")
    ejes[0].set_ylabel("Bray-Curtis dissimilarity")
    ejes[0].set_title("Same species, different strain")
    ejes[0].legend(fontsize=7)
    ejes[0].grid(alpha=0.3)
    for k, estilo in [(15, ":"), (21, "-"), (31, "--")]:
        sel = [f for f in filas if f["k"] == str(k) and f["par"] == "X-X_cepa"]
        d = [100 * float(f["divergencia"]) for f in sel]
        ejes[1].plot(d, [float(f["jaccard_simka"]) for f in sel], estilo, marker="o", color=COLOR["Simka"],
                     label="Simka, k = %d" % k)
        ejes[1].plot(d, [float(f["jaccard_nomka"]) for f in sel], estilo, marker="s", color=COLOR["Nomka"],
                     label="Nomka, k = %d" % k)
    ejes[1].set_xlabel("SNP divergence between strains X and X' (%)")
    ejes[1].set_ylabel("Jaccard dissimilarity")
    ejes[1].set_title("Same species, different strain (presence/absence)")
    ejes[1].legend(fontsize=7)
    ejes[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG + "/fig_sintetico.png", dpi=200)


def tabla_resumen(med):
    calidad = leer_csv("experimentos/tablas/calidad.csv")
    lineas = ["| Algorithm | k | N | Wall time (s) | CPU time (s) | CPU (%) | Peak RAM (MB) | Input (MB) | Temp. disk (MB) | Silhouette BC | NN same site BC (%) |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
    for k in ["21", "31"]:
        for n in ["2", "4", "6", "9"]:
            for programa in ["Simka-Python", "Nomka"]:
                filas = [f for f in med if f["programa"] == programa and f["k"] == k and f["n"] == n]
                if not filas:
                    continue
                t = media_de([float(f["tiempo_s"]) for f in filas])
                c = media_de([float(f["cpu_s"]) for f in filas])
                p = statistics.mean(float(f["cpu_pct"]) for f in filas)
                r = media_de([float(f["ram_mb"]) for f in filas])
                q = [f for f in calidad if f["k"] == k and f["n"] == n and f["distancia"] == "braycurtis"
                     and f["programa"] == ("Simka" if programa == "Simka-Python" else "Nomka")]
                sil = q[0]["silueta"] if q else "-"
                vec = q[0]["vecino_ok_pct"] if q else "-"
                lineas.append("| %s | %s | %s | %.1f ± %.1f | %.1f ± %.1f | %.0f | %.0f ± %.0f | %s | 0 | %s | %s |"
                              % (ETIQUETA[programa], k, n, t[0], t[1], c[0], c[1], p, r[0], r[1],
                                 filas[0]["entrada_mb"], sil, vec))
    with open("experimentos/tablas/resumen.md", "w") as f:
        f.write("\n".join(lineas) + "\n")
    print("\n".join(lineas))


def main():
    os.makedirs(FIG, exist_ok=True)
    med = [f for f in leer_csv("experimentos/tablas/mediciones.csv") if f["codigo_salida"] == "0"]
    figura_recursos(med)
    for distancia in ["braycurtis", "jaccard"]:
        figura_matrices(21, distancia)
    if os.path.exists("experimentos/tablas/delta_pares.csv"):
        figura_delta(21)
    if os.path.exists("experimentos/tablas/sintetico.csv"):
        figura_sintetico()
    tabla_resumen(med)


if __name__ == "__main__":
    main()
