# Nomka: comparing metagenomic samples with k-mer counts

Course project for BCOM4006 (Algorithms in Computational Biology, Universidad de
los Andes, 2026-2), based on **Simka** (Benoit et al., 2016, *PeerJ Computer
Science* 2:e94; https://github.com/GATB/simka).

Given N metagenomic samples (sets of reads), the programs compute the N × N
Bray–Curtis and Jaccard dissimilarity matrices from canonical k-mer counts,
with the same input and output formats as Simka.

| File | What it does |
|---|---|
| `src/kmer_distances.py` | In-memory algorithm: one hash table (dict) per sample, solid k-mers (count ≥ t), pairwise comparison. Gives exactly the same matrices as Simka |
| `src/Nomka.py` | Same as `kmer_distances.py` plus one step: k-mers that differ in one letter (one SNP) are grouped into a family before comparing |
| `src/Simka.py` | Reimplementation that follows the phases of the original C++ Simka (partitions, sorting count, merging count); used as an independent check |
| `datos_simka/` | Example data of Simka (samples A–E) and the official Simka output matrices for k = 21, 31 and t = 0, 2 (from `tests/truth` of GATB/simka) |
| `experimentos/` | Scripts, tables and figures of the experiments (HMP and synthetic data) |

Requirements: Python 3, standard library only. The figures also need
`matplotlib` and `scipy`.

## Usage (from the repository root)

```bash
# <samples file> <k> <t> <output prefix>; writes <prefix>_braycurtis.csv and <prefix>_jaccard.csv
python3 src/kmer_distances.py datos_simka/simka_input.txt 21 2 resultados/ejemplo
python3 src/Nomka.py          datos_simka/simka_input.txt 21 2 resultados/nomka_ejemplo
```

The samples file uses the Simka format, one sample per line:
`name: file1 , file2 ; file3` (`,` joins parts of a sample, `;` separates
paired files). FASTA and FASTQ are accepted.

## Verification on the Simka example

```bash
python3 src/Simka.py      # prints True 4 times (BC and Jaccard, t = 0 and t = 2)
cmp resultados/ejemplo_braycurtis.csv datos_simka/referencia_k21_t2/mat_abundance_braycurtis.csv
cmp resultados/ejemplo_jaccard.csv    datos_simka/referencia_k21_t2/mat_presenceAbsence_jaccard.csv
```

`cmp` prints nothing when the files are identical. Nomka also gives identical
matrices on this example, because it has no k-mers that differ in one letter.

## Experiments

Data: 9 samples of the Human Microbiome Project (SRA project SRP002163; 3 stool,
3 tongue dorsum, 3 posterior fornix), first 100,000 reads of R1. The reads are
not stored in this repository.

```bash
python3 experimentos/descargar_hmp.py 100000      # needs sratoolkit (fastq-dump); writes datos_hmp/
python3 experimentos/medir.py 3                   # time, CPU, RAM; N = 2, 4, 6, 9; 3 runs in parallel
python3 experimentos/medir.py 1 control           # control runs (each program alone)
python3 experimentos/sintetico.py                 # synthetic community with known SNP divergence
python3 experimentos/comparar_matrices.py         # quality metrics and Δ per pair
python3 experimentos/verificar.py --independiente # reproducibility, postconditions, Simka.py check
python3 -m venv .venv && .venv/bin/pip install matplotlib scipy
.venv/bin/python experimentos/graficas.py         # figures and summary table
```

Results already obtained:

- `resultados/hmp/`: output matrices of every run (`<program>_k<k>_n<N>_r<rep>_{braycurtis,jaccard}.csv`) and their logs.
- `experimentos/tablas/`:
  - `mediciones.csv`: one row per run.
  - `resumen.md`: cost and quality summary.
  - `calidad.csv`, `delta_pares.csv`, `sintetico.csv`, `verificacion.txt`.
- `experimentos/figuras/`:
  - cost vs N;
  - Bray–Curtis and Jaccard matrices with UPGMA;
  - Δ per pair;
  - synthetic experiment.

## License note

`datos_simka/` comes from the Simka repository (AGPL-3.0).
