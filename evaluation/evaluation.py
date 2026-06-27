import subprocess
import os
import random
import numpy as np
import matplotlib.pyplot as plt

def subsample_reads(input_fq, output_fq, coverage):
    """Subsample reads from input_fq to achieve the desired coverage."""
    # Read all reads from the input file
    with open(input_fq, 'r') as f:
        reads = f.read().split('\n\n')  # Split by empty lines (each read is separated by a blank line)
    
    # Calculate the number of reads needed to achieve the desired coverage
    # This is a simplified approach; real subsampling may require more sophisticated methods
    num_reads_needed = int(coverage * len(reads))
    
    # Shuffle the reads to ensure randomness
    random.shuffle(reads)
    
    # Write the subsampled reads to the output file
    with open(output_fq, 'w') as f:
        f.write('\n\n'.join(reads[:num_reads_needed]))

def calculate_n50(fasta_file):
    """Calculate the N50 contig length from a FASTA file."""
    contig_lengths = []
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        if lines[i].startswith('>'):
            i += 1
            contig_length = 0
            while i < len(lines) and not lines[i].startswith('>'):
                contig_length += len(lines[i].strip())
                i += 1
            contig_lengths.append(contig_length)
    
    contig_lengths.sort(reverse=True)
    total_length = sum(contig_lengths)
    cumulative_length = 0
    n50 = 0
    for length in contig_lengths:
        cumulative_length += length
        if cumulative_length >= total_length / 2:
            n50 = length
            break
    return n50

# Configuration parameters for Ilesta assembler
config = {
    # Read filtering and alignment parameters
    'reads_fq': 'evaluation/test_data/SRR28262566.fastq',
    'output_dir': 'evaluation/test_data/output',
    'threads': 4,
    'paf': 'alignments.paf',
    'min_read_length': 1000,
    'min_base_quality': 10.0,
    'genome_size': None,

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
coverages = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
n50_values = []

for coverage in coverages:
    # Subsample reads
    subsampled_fq = f"{config['output_dir']}/subsampled_{coverage:.1f}.fq"
    subsample_reads(config['reads_fq'], subsampled_fq, coverage)
    
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

# Plot the N50 values
coverages, n50_values = zip(*n50_values)
plt.figure(figsize=(10, 6))
plt.plot(coverages, n50_values, marker='o')
plt.xlabel('Coverage')
plt.ylabel('N50 Contig Length')
plt.title('N50 Contig Length vs Coverage')
plt.grid(True)
plt.show()