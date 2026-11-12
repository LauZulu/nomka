"""
Calculates the Bray-Curtis and Jaccard distances between N metagenomic samples,
counting their k-mers in memory with one hash table (dict) per sample.

Usage:
    python3 kmer_distances.py <samples file> <k> <t> <output prefix>

Example:
    python3 kmer_distances.py simka_input.txt 21 2 results/example

The samples file uses the same format as the input of Simka:
one sample per line, "name: file1 , file2 ; file3".
"""

import os
import sys


# ---------------------------------------------------------------------------
# 1. Reading the input
# ---------------------------------------------------------------------------

def read_samples_file(samples_file):
    """
    Reads the list of samples.
    :param samples_file: file with one sample per line ("name: file1 , file2 ; file3")
    :return: two lists with the same order: the names and the files of each sample
    """
    folder = os.path.dirname(os.path.abspath(samples_file))
    names = []
    files = []
    with open(samples_file) as f:
        for line in f:
            line = line.strip()
            if line == "":
                continue
            colon = line.index(":")
            name = line[:colon].strip()
            # Files of the same sample are separated by ',' or by ';' (paired reads)
            parts = line[colon + 1:].replace(";", ",").split(",")
            sample_files = []
            for part in parts:
                file_name = part.strip()
                if file_name != "":
                    # Relative names are relative to the folder of the samples file
                    if not os.path.isabs(file_name):
                        file_name = os.path.join(folder, file_name)
                    sample_files.append(file_name)
            names.append(name)
            files.append(sample_files)
    return names, files


def count_file(file_name, k, counts):
    """
    Reads a FASTA or FASTQ file and adds the k-mers of each read to the counts.
    Each read is counted as soon as it is read, so the file is never kept in memory.
    :param file_name: FASTA (header lines start with '>') or FASTQ (header lines start with '@')
    :param k: length of the k-mers
    :param counts: dict of the sample (canonical k-mer -> count)
    :return: number of reads in the file
    """
    num_reads = 0
    with open(file_name) as f:
        first_line = f.readline()
        if first_line.startswith("@"):
            # FASTQ: 4 lines per read (header, sequence, '+', qualities)
            line = first_line
            while line != "":
                sequence = f.readline().strip()
                f.readline()
                f.readline()
                add_read(sequence, k, counts)
                num_reads += 1
                line = f.readline()
        else:
            # FASTA: a header line and then one or more lines with the sequence
            sequence = ""
            line = first_line
            while line != "":
                if line.startswith(">"):
                    if sequence != "":
                        add_read(sequence, k, counts)
                        num_reads += 1
                    sequence = ""
                else:
                    sequence = sequence + line.strip()
                line = f.readline()
            if sequence != "":
                add_read(sequence, k, counts)
                num_reads += 1
    return num_reads


# ---------------------------------------------------------------------------
# 2. Counting the k-mers of one sample
# ---------------------------------------------------------------------------

def is_valid(kmer):
    """
    :return: True if the k-mer only has the letters A, C, G and T (no N)
    """
    for c in kmer:
        if c != "A" and c != "C" and c != "G" and c != "T":
            return False
    return True


def reverse_complement(kmer):
    """
    :return: reverse complement of a DNA sequence (for example, AAC returns GTT)
    """
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    answer = ""
    for i in range(len(kmer) - 1, -1, -1):
        answer = answer + complement[kmer[i]]
    return answer


def canonical(kmer):
    """
    :return: the alphabetically smaller between the k-mer and its reverse complement
    """
    rc = reverse_complement(kmer)
    if kmer <= rc:
        return kmer
    return rc


def add_read(read, k, counts):
    """
    Adds the canonical k-mers of one read to the hash table of the sample.
    :param read: sequence of the read
    :param k: length of the k-mers
    :param counts: dict where key = canonical k-mer, value = number of times it appears
    """
    read = read.upper()
    for p in range(len(read) - k + 1):
        kmer = read[p:p + k]
        if is_valid(kmer):
            c = canonical(kmer)
            counts[c] = counts.get(c, 0) + 1


def remove_rare_kmers(counts, t):
    """
    Removes the k-mers that appear less than t times (keeps the solid k-mers).
    :param counts: k-mer counts of one sample
    :param t: minimum abundance
    """
    to_remove = []
    for kmer in counts:
        if counts[kmer] < t:
            to_remove.append(kmer)
    for kmer in to_remove:
        del counts[kmer]


# ---------------------------------------------------------------------------
# 3. Comparing the samples
# ---------------------------------------------------------------------------

def distance_matrices(tables):
    """
    Calculates the Bray-Curtis and Jaccard distance matrices.
    :param tables: tables[i] has the solid k-mer counts of sample i
    :return: the two N x N matrices (lists of lists)
    """
    n = len(tables)
    bray_curtis = [[0.0] * n for _ in range(n)]
    jaccard = [[0.0] * n for _ in range(n)]

    # Total number of solid k-mers (C) and number of distinct solid k-mers (d) of each sample
    total = []
    distinct = []
    for i in range(n):
        total.append(sum(tables[i].values()))
        distinct.append(len(tables[i]))

    for i in range(n):
        for j in range(i + 1, n):
            # The smaller table is traversed and each k-mer is searched in the larger one
            small = tables[i]
            big = tables[j]
            if len(small) > len(big):
                small = tables[j]
                big = tables[i]
            shared = 0   # sum of min(count in i, count in j): numerator of Bray-Curtis
            common = 0   # number of distinct k-mers present in both: numerator of Jaccard
            for kmer in small:
                if kmer in big:
                    shared += min(small[kmer], big[kmer])
                    common += 1
            bc = 1.0 - 2.0 * shared / (total[i] + total[j])
            jac = 1.0 - common / (distinct[i] + distinct[j] - common)
            bray_curtis[i][j] = bc
            bray_curtis[j][i] = bc
            jaccard[i][j] = jac
            jaccard[j][i] = jac
    return bray_curtis, jaccard


def write_matrix(file_name, names, matrix):
    """
    Writes a distance matrix with the same format as Simka (separator ';', 6 decimals).
    """
    with open(file_name, "w") as out:
        out.write(";" + ";".join(names) + "\n")
        for i in range(len(names)):
            values = ["%.6f" % matrix[i][j] for j in range(len(names))]
            out.write(names[i] + ";" + ";".join(values) + "\n")


# ---------------------------------------------------------------------------
# Main program
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 5:
        print("Usage: python3 kmer_distances.py <samples file> <k> <t> <output prefix>")
        sys.exit(1)
    samples_file = sys.argv[1]
    k = int(sys.argv[2])
    t = int(sys.argv[3])
    output_prefix = sys.argv[4]

    # 1. Read the list of samples
    names, files = read_samples_file(samples_file)

    # 2. Count the k-mers of each sample and keep only the solid k-mers
    tables = []
    print("Sample\tReads\tDistinctSolidKmers\tTotalSolidKmers")
    for i in range(len(names)):
        counts = {}
        num_reads = 0
        for file_name in files[i]:
            num_reads += count_file(file_name, k, counts)
        remove_rare_kmers(counts, t)
        tables.append(counts)
        print(names[i] + "\t" + str(num_reads) + "\t" + str(len(counts)) + "\t" + str(sum(counts.values())))

    # 3. Calculate the distance matrices and write them
    bray_curtis, jaccard = distance_matrices(tables)
    write_matrix(output_prefix + "_braycurtis.csv", names, bray_curtis)
    write_matrix(output_prefix + "_jaccard.csv", names, jaccard)
    print("Matrices written with prefix " + output_prefix)


if __name__ == "__main__":
    main()
