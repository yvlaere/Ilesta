import subprocess
import os
import random
import numpy as np
import matplotlib.pyplot as plt

def subsample_reads(input_fq, output_fq, coverage, genome_size=6741160, seed=None):
    """
    Randomly subsample a FASTQ file to approximately the requested coverage.

    Parameters
    ----------
    input_fq : str
        Input FASTQ file.
    output_fq : str
        Output FASTQ file.
    coverage : float
        Desired sequencing coverage (e.g. 5, 10, 30).
    genome_size : int
        Genome size in base pairs.
    seed : int, optional
        Random seed for reproducibility.
    """

    rng = random.Random(seed)

    target_bp = genome_size * coverage

    # ------------------------------------------------------------------
    # First pass: determine the total number of bases in the FASTQ
    # ------------------------------------------------------------------
    total_bp = 0

    with open(input_fq) as fin:
        while True:
            header = fin.readline()
            if not header:
                break

            seq = fin.readline()
            sep = fin.readline()
            qual = fin.readline()

            if not (seq and sep and qual):
                raise ValueError("Malformed FASTQ: incomplete record.")

            total_bp += len(seq.rstrip())

    if total_bp == 0:
        raise ValueError("Input FASTQ contains no reads.")

    keep_probability = min(1.0, target_bp / total_bp)

    # ------------------------------------------------------------------
    # Second pass: randomly retain reads
    # ------------------------------------------------------------------
    with open(input_fq) as fin, open(output_fq, "w") as fout:
        while True:
            header = fin.readline()
            if not header:
                break

            seq = fin.readline()
            sep = fin.readline()
            qual = fin.readline()

            if not (seq and sep and qual):
                raise ValueError("Malformed FASTQ: incomplete record.")

            if rng.random() < keep_probability:
                fout.write(header)
                fout.write(seq)
                fout.write(sep)
                fout.write(qual)

def get_contig_lengths(fasta_file):
    """Return a list containing the length of every contig in a FASTA file."""
    contig_lengths = []

    with open(fasta_file, "r") as f:
        contig_length = 0
        for line in f:
            if line.startswith(">"):
                if contig_length:
                    contig_lengths.append(contig_length)
                contig_length = 0
            else:
                contig_length += len(line.strip())

        if contig_length:
            contig_lengths.append(contig_length)

    return contig_lengths

def calculate_largest_contig(fasta_file):
    contig_lengths = get_contig_lengths(fasta_file)
    return max(contig_lengths) if contig_lengths else 0

def calculate_n50(fasta_file):
    contig_lengths = sorted(get_contig_lengths(fasta_file), reverse=True)

    total_length = sum(contig_lengths)
    cumulative = 0

    for length in contig_lengths:
        cumulative += length
        if cumulative >= total_length / 2:
            return length

    return 0

# Configuration parameters for Ilesta assembler
config = {
    # Read filtering and alignment parameters
    'reads_fq': 'evaluation/test_data/SRR28262566.fastq',
    'output_dir': 'evaluation/output',
    'threads': 4,
    'paf': 'alignments.paf',
    'min_read_length': 1000,
    'min_base_quality': 10.0,
    'genome_size': 6741160,

    # Alignment filtering parameters
    'min_overlap_length': 2000,
    'min_overlap_count': 3,
    'min_percent_identity': 5.0,
    'overhang_ratio': 0.8,

    # Assembly parameters
    'output_prefix': 'unitigs',
    'max_bubble_length': 100,
    'min_support_ratio': 1.1,
    'max_tip_len': 4,
    'fuzz': 10,
    'cleanup_iterations': 3,
    'short_edge_ratio': 0.8,
    'completion_enabled': False,
    'completion_rounds': 1,
    'completion_min_alignment_len': 2000,
    'completion_min_identity': 0.8,
}

# Path to the Ilesta binary
ilesta_path = 'target/debug/Ilesta'

# Create output directory if it doesn't exist
os.makedirs(config['output_dir'], exist_ok=True)

# Define coverages for subsampling
coverages = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
n50_values = []
largest_contig_values = []

for coverage in coverages:
    # Subsample reads
    subsampled_fq = f"{config['output_dir']}/subsampled_{coverage:.1f}.fq"
    subsample_reads(config['reads_fq'], subsampled_fq, coverage)
    # print the number of reads in the subsampled file
    with open(subsampled_fq, 'r') as f:
        num_reads = sum(1 for line in f) // 4
    print(f"Subsampled {num_reads} reads for coverage {coverage:.1f}x")
    
    # Run the Ilesta assembler for the subsampled reads
    output_dir = f"{config['output_dir']}/coverage_{coverage:.1f}"
    os.makedirs(output_dir, exist_ok=True)
    
    # Modify the configuration for the current subsampled set
    config['output_dir'] = output_dir
    config['reads_fq'] = subsampled_fq
    
    # Run the Ilesta assembler
    subprocess.run([ilesta_path, 'assemble',
                    f'--reads-fq={config["reads_fq"]}',
                    f'--output-dir={config["output_dir"]}',
                    f'--threads={config["threads"]}',
                    f'--paf={config["paf"]}',
                    f'--min-read-length={config["min_read_length"]}',
                    f'--min-base-quality={config["min_base_quality"]}',
                    f'--genome-size={config["genome_size"]}',
                    f'--min-overlap-length={config["min_overlap_length"]}',
                    f'--min-overlap-count={config["min_overlap_count"]}',
                    f'--min-percent-identity={config["min_percent_identity"]}',
                    f'--overhang-ratio={config["overhang_ratio"]}',
                    f'--output-prefix={config["output_prefix"]}',
                    f'--max-bubble-length={config["max_bubble_length"]}',
                    f'--min-support-ratio={config["min_support_ratio"]}',
                    f'--max-tip-len={config["max_tip_len"]}',
                    f'--fuzz={config["fuzz"]}',
                    f'--cleanup-iterations={config["cleanup_iterations"]}',
                    f'--short-edge-ratio={config["short_edge_ratio"]}',
                    '--completion-enabled' if config["completion_enabled"] is False else '--completion-enabled=true',
                    f'--completion-rounds={config["completion_rounds"]}',
                    f'--completion-min-alignment-len={config["completion_min_alignment_len"]}',
                    f'--completion-min-identity={config["completion_min_identity"]}'],
                   check=True)
    
    # Calculate N50 for the current assembly
    n50 = calculate_n50(f"{output_dir}/{config['output_prefix']}.fa")
    n50_values.append((coverage, n50))

    # Calculate the largest contig for the current assembly
    largest_contig = calculate_largest_contig(f"{output_dir}/{config['output_prefix']}.fa")
    largest_contig_values.append((coverage, largest_contig))

# Plot the N50 values
coverages, n50_values = zip(*n50_values)
plt.figure(figsize=(10, 6))
plt.plot(coverages, n50_values, marker='o')
plt.xlabel('Coverage')
plt.ylabel('N50 Contig Length')
plt.title('N50 Contig Length vs Coverage')
plt.grid(True)
plt.tight_layout()
plt.savefig("evaluation/output/n50_vs_coverage.png", dpi=300)
plt.close()

# Plot the largest contig values
coverages, largest_contig_values = zip(*largest_contig_values)
plt.figure(figsize=(10, 6))
plt.plot(coverages, largest_contig_values, marker='o', color='orange')
plt.xlabel('Coverage')
plt.ylabel('Largest Contig Length')
plt.title('Largest Contig Length vs Coverage')
plt.grid(True)
plt.tight_layout()
plt.savefig("evaluation/output/largest_contig_vs_coverage.png", dpi=300)