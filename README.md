# Nomka

Proyecto de BCOM4006 (Algoritmos en Biología Computacional, 2026-2) sobre Simka
(Benoit et al., 2016). Se calculan distancias Bray-Curtis y Jaccard entre muestras
metagenómicas a partir de k-mers.

- `src/kmer_distances.py`: una tabla hash por muestra; da lo mismo que Simka.
- `src/Nomka.py`: igual, pero junta los k-mers que difieren en una letra (SNP).
- `src/Simka.py`: versión en Python de las fases de Simka, para comprobar.
- `experimentos/`: scripts, tablas y figuras de las pruebas con HMP.
- `resultados/hmp/`: matrices de cada corrida.
- `datos_simka/`: ejemplo de Simka con sus matrices de referencia.

## Cómo correrlo

Solo necesita Python 3.

```bash
python3 src/kmer_distances.py datos_simka/simka_input.txt 21 2 resultados/ejemplo
python3 src/Nomka.py datos_simka/simka_input.txt 21 2 resultados/nomka_ejemplo
```

Los argumentos son: archivo de muestras, k, t y prefijo de salida.
