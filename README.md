# Ilesta

[![DOI](https://zenodo.org/badge/1043754381.svg)](https://doi.org/10.5281/zenodo.18699305)

Ilesta assembles long reads (ONT or PacBio) using the Overlap-Layout-Consensus approach. It finds overlaps between reads, builds and simplifies an overlap graph, and calls a consensus sequence for the final contigs.

[Install](#installation) · [Quick usage](#quick-usage) · [How it works](#how-it-works) · [Command line usage](#command-line-usage)

## Installation

Ilesta uses [minimap2](https://github.com/lh3/minimap2) for read alignment. The Bioconda package installs minimap2 automatically. For binary and source installations, install minimap2 separately and make sure it is available on your `PATH`.

### Option 1: Bioconda (recommended)

#### Using Conda

```bash
conda create -n ilesta_env -c conda-forge -c bioconda ilesta
conda activate ilesta_env
```

#### Using Pixi

```bash
mkdir ilesta_workspace && cd ilesta_workspace
pixi init
pixi add -c conda-forge -c bioconda ilesta
```

### Option 2: Precompiled binary

```bash
mkdir Ilesta && cd Ilesta
wget https://github.com/yvlaere/Ilesta/releases/latest/download/ilesta-linux-x86_64
chmod +x ilesta-linux-x86_64
mv ilesta-linux-x86_64 Ilesta

# Optional: add Ilesta to PATH
export PATH="$PWD:$PATH"
```

### Option 3: Build from source

Use a Rust toolchain that satisfies the requirements in `Cargo.toml`.

```bash
git clone https://github.com/yvlaere/Ilesta.git
cd Ilesta
cargo build --release

# Optional: add the binary to PATH
export PATH="$PWD/target/release:$PATH"
```

## Quick usage

> [!CAUTION]
> For Ilesta versions earlier than 1.5.0, use `--overhang-ratio 0.2` when running `assemble`. This improves assembly completeness, especially in datasets with lower coverage or shorter reads.

```bash
Ilesta assemble \
  --reads-fq reads.fastq \
  --output-dir out_dir \
  --threads 16 \
  --min-read-length 1000 \
  --min-base-quality 10
```

The read-length and mean-quality filters are optional and read filtering is **disabled by default**. FASTQ input may be plain or gzip-compressed. Ilesta estimates genome size from read overlaps unless `--genome-size` is supplied. Use `--read-type` to select `ont` (default), `pb-clr` or `pb-hifi`.

### Output files

| File | Description |
| --- | --- |
| `contigs.fa` | Final assembled contigs in FASTA format. |
| `contigs.gfa` | Final assembly graph in GFA format. |

Intermediate files (alignments, selected reads, overlaps) are removed by default. Add `--keep-intermediate` to keep them.

## How it works

1. **Prepare and align reads.** Optional length and quality filters are applied, then reads are selected up to the target coverage. Ilesta estimates genome size from the overlaps unless you provide it, and uses minimap2 to find read-to-read overlaps.
2. **Build the overlap graph.** Filtered overlaps become an oriented graph of reads. Ilesta simplifies the graph by removing transitive and redundant edges, bubbles, and tips, then compresses graph paths into draft contigs.
3. **Refine contigs.** By default, Ilesta tries to rescue reads that are poorly represented in the draft. It removes contained contigs, uses read evidence to bridge contigs, and detects circular contigs from overlapping ends.
4. **Polish and write output.** By default, reads are mapped back to contigs for windowed consensus polishing. The final assembly is written as `contigs.fa` and `contigs.gfa`.

## Command Line Usage

### `Ilesta --help`

```bash
Ilesta --help
```

```text
De novo genome assembly for long reads using an overlap graph

Usage: Ilesta <COMMAND>

Commands:
  align                Read selection and alignment
  alignment-filtering  Alignment filtering
  assemble             Full genome assembly pipeline
  scaffold             Read-guided bridging of an existing contig FASTA
  help                 Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```

### `Ilesta assemble --help`

```bash
Ilesta assemble --help
```

```text
Full genome assembly pipeline

Usage: Ilesta assemble [OPTIONS] --reads-fq <READS_FQ>

Options:
  -h, --help  Print help

Output:
  -p, --output-prefix <OUTPUT_PREFIX> Output prefix [default: contigs]
  -o, --output-dir <OUTPUT_DIR>        Output directory [default: .]
      --write-dot                      Also write the overlap graph as graph.dot
      --keep-intermediate              Keep intermediate files (alignments, selected reads, overlaps) in the output directory instead of removing them; by default only the FASTA and GFA remain

Read processing and alignment:
  -r, --reads-fq <READS_FQ>
          Input reads in FASTQ format
      --read-type <READ_TYPE>
          Read type / sequencing technology [default: ont] [possible values: ont, pb-clr, pb-hifi]
  -t, --threads <THREADS>
          Number of threads [default: 4]
  -a, --paf <PAF>
          Output PAF filename [default: alignments.paf]
      --min-read-length <MIN_READ_LENGTH>
          Minimum read length (0 = no filtering) [default: 0]
  -q, --min-base-quality <MIN_BASE_QUALITY>
          Minimum average quality (0 = no filtering) [default: 0]
      --genome-size <GENOME_SIZE>
          Optional input genome size (if not provided, will be estimated from data)
      --target-coverage <TARGET_COVERAGE>
          Target read coverage for assembly input downsampling [default: 50]
      --minimap2-k <MINIMAP_BATCH_SIZE>
          Override minimap2's -K query minibatch size for self-alignment (e.g. 500M) [default: 500M]
      --minimap2-f <MINIMAP2_F>
          Override minimap2's -f high-frequency minimizer filter
      --minimap2-u <MINIMAP2_U>
          Override minimap2's -U minimizer occurrence bounds [default: 10,5000]

Alignment filtering:
  -l, --min-overlap-length <MIN_OVERLAP_LENGTH>
          Minimum overlap length [default: 2000]
  -c, --min-overlap-count <MIN_OVERLAP_COUNT>
          Minimum overlap count [default: 3]
  -i, --min-percent-identity <MIN_PERCENT_IDENTITY>
          Minimum percent identity [default: 5]
      --overhang-ratio <OVERHANG_RATIO>
          Overhang ratio [default: 0.2]
      --overlaps <OVERLAPS>
          Pre-computed overlaps binary file (optional, if provided skips alignment filtering)

Assembly:
      --max-bubble-length <MAX_BUBBLE_LENGTH>
          Maximum bubble length (used during bubble removal) [default: 100]
      --min-support-ratio <MIN_SUPPORT_RATIO>
          Minimum support ratio for bubble removal [default: 1.1]
      --max-tip-len <MAX_TIP_LEN>
          Maximum tip length for tip trimming [default: 4]
      --fuzz <FUZZ>
          Fuzz parameter for transitive edge reduction [default: 10]
      --cleanup-iterations <CLEANUP_ITERATIONS>
          Number of cleanup iterations to run [default: 3]
      --short-edge-ratio <SHORT_EDGE_RATIO>
          Short edge removal ratio (heuristic simplification) [default: 0.8]
      --rescue-min-overlap-length <RESCUE_MIN_OVERLAP_LENGTH>
          Minimum overlap length in the rescue stage (the unmapped reads are short and low coverage) [default: 2000]
      --no-rescue
          Disable rescue of poorly represented reads (enabled by default)
      --seed <SEED>
          Seed (the same seed gives the same result) [default: 42]

Scaffolding:
      --scaffold-rounds <ROUNDS>
          Contig bridging rounds (read-guided scaffolding of the final contigs; 0 disables) [default: 2]
      --scaffold-min-reads <MIN_READS>
          Minimum number of reads spanning a bridge between two contigs [default: 3]
      --scaffold-dominance <DOMINANCE>
          A competing bridge at the same contig end may have at most this fraction of the best bridge's reads [default: 0.3]
      --scaffold-min-alignment-len <MIN_ALIGNMENT_LEN>
          Minimum read-to-contig alignment length used for bridging [default: 1000]
      --scaffold-min-identity <MIN_IDENTITY>
          Minimum read-to-contig alignment identity used for bridging [default: 0.8]

Consensus:
      --consensus-rounds <CONSENSUS_ROUNDS>
          Consensus polishing rounds on the final contigs (0 disables) [default: 2]
      --consensus-window-len <CONSENSUS_WINDOW_LEN>
          Window length (bp) for consensus calling [default: 500]
      --consensus-min-depth <CONSENSUS_MIN_DEPTH>
          Minimum read segments per window to re-call it (else keep draft) [default: 3]
      --consensus-max-depth <CONSENSUS_MAX_DEPTH>
          Maximum read segments used per window [default: 25]
      --consensus-repair
          Replace stretches that no read supports but that many reads bridge (aggressive polishing)
```

## Development

Run the default test suite with:

```bash
cargo test
```