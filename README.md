# Ilesta

[![DOI](https://zenodo.org/badge/1043754381.svg)](https://doi.org/10.5281/zenodo.18699305)

Ilesta is a de novo genome assembler for long reads. It processes all-vs-all alignments of long-read sequencing data to detect overlaps, construct an overlap graph, and generate assembly unitigs. It does not have a consensus step, so it should be paired with a polishing tool.

## Installation

Option 1: Bioconda (recommended)

```bash
conda create -n Ilesta_env
conda activate Ilesta_env
conda install -c bioconda ilesta
```

Option 2: Precompiled binaries:

```bash
# Download the binary
mkdir Ilesta; cd Ilesta
wget https://github.com/yvlaere/Ilesta/releases/download/v1.0.0/ilesta-linux-x86_64
chmod +x ilesta-linux-x86_64
mv ilesta-linux-x86_64 Ilesta

# Add Ilesta to your PATH (optional, for easy access):
export PATH="$PWD:$PATH"
```

Option 3: Build from source (requires Rust 1.85.0 or above):

```bash
# Build the release
git clone https://github.com/yvlaere/Ilesta.git
cd Ilesta
cargo build --release

# Add Ilesta to your PATH (optional, for easy access):
export PATH="$PWD/target/release:$PATH"
```

## Usage

A demonstration of the usage of Ilesta using data from: https://doi.org/10.3389/fmicb.2025.1532788

```bash
# Download data, using SRR28262566 as an example
prefetch SRR28262566; cd SRR28262566
fasterq-dump SRR28262566.sra

# Filter reads
seqkit seq -m 10000 -Q 10 SRR28262566.fastq -o filtered.fq

# All-to-all alignment of the reads
minimap2 -x ava-ont -t4 filtered.fq filtered.fq > overlapped_reads.paf

# Start the assembly
Ilesta assemble --input-paf overlapped_reads.paf --reads-fq filtered.fq -o out_dir
```

This will produce:
- `out_dir/unitigs.fa` (unitigs in FASTA format)
- `out_dir/unitigs.gfa` (assembly graph in GFA format)
- `out_dir/graph.dot` (overlap graph visualization)

```bash
# Long read polishing
minipolish filtered.fq out_dir/unitigs.gfa > polished.gfa

# visualize the assembly graph
Bandage image out_dir/unitigs.gfa out_dir/unitigs.png
```
Ilesta sucesfully assembles the whole bacterial chromosome as one unitig. The smaller unitigs are plasmids. 
![Bandage visualization](image.png)

```bash
# visualize the overlap graph (not recommended, the graph is usually too large for convenient visualization)
dot -Tpng out_dir/graph.dot -o out_dir/graph.png
```

## Command Line Usage

```
Ilesta --help

Usage: Ilesta <COMMAND>

Commands:
  alignment-filtering  Alignment filtering
  assemble             Full genome assembly pipeline
  help                 Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```
```
Ilesta assemble --help

Usage: Ilesta assemble [OPTIONS] --reads-fq <READS_FQ>

Options:
  -f, --input-paf <INPUT_PAF>
          Input PAF file (optional if --overlaps is provided)
  -l, --min-overlap-length <MIN_OVERLAP_LENGTH>
          Minimum overlap length [default: 2000]
  -c, --min-overlap-count <MIN_OVERLAP_COUNT>
          Minimum overlap count [default: 3]
  -i, --min-percent-identity <MIN_PERCENT_IDENTITY>
          Minimum percent identity [default: 5]
      --overhang-ratio <OVERHANG_RATIO>
          Overhang ratio [default: 0.8]
      --overlaps <OVERLAPS>
          Pre-computed overlaps binary file (optional, if provided skips alignment filtering)
  -r, --reads-fq <READS_FQ>
          Input reads in FASTQ format
  -p, --output-prefix <OUTPUT_PREFIX>
          Output prefix [default: unitigs]
  -o, --output-dir <OUTPUT_DIR>
          Output directory [default: .]
      --max-bubble-length <MAX_BUBBLE_LENGTH>
          Maximum bubble length (used during bubble removal) [default: 100]
      --min-support-ratio <MIN_SUPPORT_RATIO>
          Minimum support ratio for bubble removal [default: 1.1]
      --max-tip-len <MAX_TIP_LEN>
          Maximum tip length for tip trimming [default: 4]
      --fuzz <FUZZ>
          Fuzz parameter for transitive edge reduction [default: 10]
      --cleanup-iterations <CLEANUP_ITERATIONS>
          Number of cleanup iterations to run [default: 2]
      --short-edge-ratio <SHORT_EDGE_RATIO>
          Short edge removal ratio (heuristic simplification) [default: 0.8]
  -h, --help
          Print help
```

## Development

Ilesta is under active development.
