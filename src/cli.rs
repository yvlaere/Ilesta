use clap::{Args, Parser, Subcommand, ValueEnum};

#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum ReadType {
    #[default]
    #[value(name = "ont")]
    Ont,
    #[value(name = "pb-clr")]
    PbClr,
    #[value(name = "pb-hifi")]
    PbHifi,
}

impl ReadType {
    pub fn minimap2_args(&self) -> &'static [&'static str] {
        match self {
            ReadType::Ont => &["-x", "ava-ont"],
            ReadType::PbClr => &["-x", "ava-pb"],
            ReadType::PbHifi => &[
                "-x", "ava-ont", "-k", "21", "-w", "11", "-g", "1000", "-m", "200", "-r", "2000",
            ],
        }
    }

    pub fn mapping_args(&self) -> &'static [&'static str] {
        match self {
            ReadType::Ont => &["-x", "map-ont"],
            ReadType::PbClr => &["-x", "map-pb"],
            ReadType::PbHifi => &["-x", "map-hifi"],
        }
    }
}

#[derive(Parser)]
#[command(
    name = "Ilesta",
    version,
    about = "De novo genome assembly for long reads using an overlap graph"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Read filtering and alignment
    Align(AlignReadsArgs),

    /// Alignment filtering
    AlignmentFiltering(AlignmentFilteringArgs),

    /// Full genome assembly pipeline
    Assemble(AssembleArgs),

    /// Read-guided bridging of an existing contig FASTA
    Scaffold(ScaffoldArgs),
}

/// Read filtering and all-vs-all alignment options (shared by `align` and `assemble`).
#[derive(Args, Clone)]
pub struct ReadOpts {
    /// Input reads in FASTQ format
    #[arg(short = 'r', long, help_heading = "Read filtering and alignment")]
    pub reads_fq: String,

    /// Read type / sequencing technology
    #[arg(long, default_value = "ont", help_heading = "Read filtering and alignment")]
    pub read_type: ReadType,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = 4, help_heading = "Read filtering and alignment")]
    pub threads: usize,

    /// Output PAF filename
    #[arg(short = 'a', long, default_value = "alignments.paf", help_heading = "Read filtering and alignment")]
    pub paf: String,

    /// Minimum read length
    #[arg(long, default_value_t = 1000, help_heading = "Read filtering and alignment")]
    pub min_read_length: u32,

    /// Minimum average quality
    #[arg(short = 'q', long, default_value_t = 10.0, help_heading = "Read filtering and alignment")]
    pub min_base_quality: f32,

    /// Optional input genome size (if not provided, will be estimated from data)
    #[arg(long, help_heading = "Read filtering and alignment")]
    pub genome_size: Option<u32>,

    /// Target read coverage for assembly input downsampling
    #[arg(long, default_value_t = 50u32, help_heading = "Read filtering and alignment")]
    pub target_coverage: u32,

    /// Minimap2 query minibatch size for self-alignment (e.g. 500M)
    #[arg(long, default_value = "500M", help_heading = "Read filtering and alignment")]
    pub minimap_batch_size: String,

    /// Override minimap2's -f high-frequency minimizer filter
    #[arg(long, help_heading = "Read filtering and alignment")]
    pub minimap2_f: Option<String>,

    /// Override minimap2's -U minimizer occurrence bounds
    #[arg(long, help_heading = "Read filtering and alignment")]
    pub minimap2_u: Option<String>,
}

/// Overlap filtering options (shared by `alignment-filtering` and `assemble`).
#[derive(Args, Clone)]
pub struct OverlapOpts {
    /// Minimum overlap length
    #[arg(short = 'l', long, default_value_t = 2000, help_heading = "Alignment filtering")]
    pub min_overlap_length: u32,

    /// Minimum overlap count
    #[arg(short = 'c', long, default_value_t = 3, help_heading = "Alignment filtering")]
    pub min_overlap_count: u32,

    /// Minimum percent identity
    #[arg(short = 'i', long, default_value_t = 5.0, help_heading = "Alignment filtering")]
    pub min_percent_identity: f32,

    /// Overhang ratio
    #[arg(long, default_value_t = 0.8, help_heading = "Alignment filtering")]
    pub overhang_ratio: f32,
}

#[derive(Args)]
pub struct AlignReadsArgs {
    /// Output directory
    #[arg(short = 'o', long, default_value = ".")]
    pub output_dir: String,

    #[command(flatten)]
    pub reads: ReadOpts,
}

#[derive(Args)]
pub struct AlignmentFilteringArgs {
    /// Output directory
    #[arg(short = 'o', long, default_value = ".")]
    pub output_dir: String,

    /// Input PAF file
    #[arg(short = 'f', long)]
    pub paf: String,

    /// Output overlaps binary file
    #[arg(long, default_value = "overlaps.bin")]
    pub output_overlaps: String,

    #[command(flatten)]
    pub overlap: OverlapOpts,

    /// Optional random seed for reproducible results
    #[arg(long)]
    pub seed: Option<u64>,
}

#[derive(Args)]
pub struct AssembleArgs {
    /// Output prefix
    #[arg(short = 'p', long, default_value = "unitigs", help_heading = "Output")]
    pub output_prefix: String,

    /// Output directory
    #[arg(short = 'o', long, default_value = ".", help_heading = "Output")]
    pub output_dir: String,

    #[command(flatten)]
    pub reads: ReadOpts,

    #[command(flatten)]
    pub overlap: OverlapOpts,

    /// Pre-computed overlaps binary file (optional, if provided skips alignment filtering)
    #[arg(long, help_heading = "Alignment filtering")]
    pub overlaps: Option<String>,

    /// Maximum bubble length (used during bubble removal)
    #[arg(long, default_value_t = 100u32, help_heading = "Assembly")]
    pub max_bubble_length: u32,

    /// Minimum support ratio for bubble removal
    #[arg(long, default_value_t = 1.1f64, help_heading = "Assembly")]
    pub min_support_ratio: f64,

    /// Maximum tip length for tip trimming
    #[arg(long, default_value_t = 4u32, help_heading = "Assembly")]
    pub max_tip_len: u32,

    /// Fuzz parameter for transitive edge reduction
    #[arg(long, default_value_t = 10u32, help_heading = "Assembly")]
    pub fuzz: u32,

    /// Number of cleanup iterations to run
    #[arg(long, default_value_t = 3u32, help_heading = "Assembly")]
    pub cleanup_iterations: u32,

    /// Short edge removal ratio (heuristic simplification)
    #[arg(long, default_value_t = 0.8f64, help_heading = "Assembly")]
    pub short_edge_ratio: f64,

    /// Disable rescue of poorly represented reads (enabled by default)
    #[arg(long, action = clap::ArgAction::SetTrue, help_heading = "Assembly")]
    pub no_rescue: bool,

    /// Also write the overlap graph as graph.dot
    #[arg(long, action = clap::ArgAction::SetTrue, help_heading = "Output")]
    pub write_dot: bool,

    #[command(flatten)]
    pub scaffold: ScaffoldOpts,

    /// Consensus polishing rounds on the final unitigs (0 disables)
    #[arg(long, default_value_t = 2u32, help_heading = "Consensus")]
    pub consensus_rounds: u32,

    /// Window length (bp) for consensus calling
    #[arg(long, default_value_t = 500usize, help_heading = "Consensus")]
    pub consensus_window_len: usize,

    /// Minimum read segments per window to re-call it (else keep draft)
    #[arg(long, default_value_t = 3usize, help_heading = "Consensus")]
    pub consensus_min_depth: usize,

    /// Maximum read segments used per window
    #[arg(long, default_value_t = 25usize, help_heading = "Consensus")]
    pub consensus_max_depth: usize,

    /// Optional random seed for reproducible results
    #[arg(long, help_heading = "Assembly")]
    pub seed: Option<u64>,
}

impl AssembleArgs {
    pub fn rescue_enabled(&self) -> bool {
        !self.no_rescue
    }
}

/// Contig bridging options (shared by `assemble` and `scaffold`).
#[derive(Args, Clone)]
pub struct ScaffoldOpts {
    /// Contig bridging rounds (read-guided scaffolding of the final unitigs; 0 disables)
    #[arg(long = "scaffold-rounds", default_value_t = 2u32, help_heading = "Scaffolding")]
    pub rounds: u32,

    /// Minimum number of reads spanning a bridge between two contigs
    #[arg(long = "scaffold-min-reads", default_value_t = 3usize, help_heading = "Scaffolding")]
    pub min_reads: usize,

    /// A competing bridge at the same contig end may have at most this fraction of the best bridge's reads
    #[arg(long = "scaffold-dominance", default_value_t = 0.3f64, help_heading = "Scaffolding")]
    pub dominance: f64,

    /// Minimum read-to-contig alignment length used for bridging
    #[arg(long = "scaffold-min-alignment-len", default_value_t = 1000u32, help_heading = "Scaffolding")]
    pub min_alignment_len: u32,

    /// Minimum read-to-contig alignment identity used for bridging
    #[arg(long = "scaffold-min-identity", default_value_t = 0.8f64, help_heading = "Scaffolding")]
    pub min_identity: f64,
}

#[derive(Args)]
pub struct ScaffoldArgs {
    /// Contigs in FASTA format (`topology=c` in the header marks circular contigs)
    #[arg(short = 'c', long)]
    pub contigs: String,

    /// Reads in FASTQ format
    #[arg(short = 'r', long)]
    pub reads_fq: String,

    /// Output FASTA
    #[arg(short = 'o', long)]
    pub output: String,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = 4usize)]
    pub threads: usize,

    /// Read type / sequencing technology
    #[arg(long, default_value = "ont")]
    pub read_type: ReadType,

    /// Consensus polishing rounds after bridging (0 = off; contig names are replaced by unitig_N)
    #[arg(long, default_value_t = 0u32)]
    pub polish_rounds: u32,

    #[command(flatten)]
    pub scaffold: ScaffoldOpts,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(args: &[&str]) -> Result<Cli, clap::Error> {
        Cli::try_parse_from(["Ilesta"].iter().chain(args))
    }

    fn assemble(args: &[&str]) -> AssembleArgs {
        match parse(args).unwrap().command {
            Commands::Assemble(a) => a,
            _ => panic!("expected assemble"),
        }
    }

    fn align(args: &[&str]) -> AlignReadsArgs {
        match parse(args).unwrap().command {
            Commands::Align(a) => a,
            _ => panic!("expected align"),
        }
    }

    #[test]
    fn minimap2_args_for_read_types() {
        assert_eq!(ReadType::Ont.minimap2_args(), &["-x", "ava-ont"]);
        assert_eq!(ReadType::PbClr.minimap2_args(), &["-x", "ava-pb"]);
        assert_eq!(ReadType::PbHifi.minimap2_args()[..2], ["-x", "ava-ont"]);
        assert_eq!(ReadType::Ont.mapping_args(), &["-x", "map-ont"]);
        assert_eq!(ReadType::PbClr.mapping_args(), &["-x", "map-pb"]);
        assert_eq!(ReadType::PbHifi.mapping_args(), &["-x", "map-hifi"]);
    }

    #[test]
    fn read_type_option() {
        assert_eq!(align(&["align", "-r", "reads.fq"]).reads.read_type, ReadType::Ont);
        assert_eq!(assemble(&["assemble", "-r", "reads.fq"]).reads.read_type, ReadType::Ont);
        for (flag, expected) in [("ont", ReadType::Ont), ("pb-clr", ReadType::PbClr), ("pb-hifi", ReadType::PbHifi)] {
            assert_eq!(align(&["align", "-r", "r.fq", "--read-type", flag]).reads.read_type, expected);
            assert_eq!(assemble(&["assemble", "-r", "r.fq", "--read-type", flag]).reads.read_type, expected);
        }
        assert!(parse(&["align", "-r", "reads.fq", "--read-type", "illumina"]).is_err());
    }

    #[test]
    fn read_options_defaults_and_overrides() {
        let a = align(&["align", "-r", "reads.fq"]);
        assert_eq!((a.reads.target_coverage, a.reads.minimap_batch_size.as_str()), (50, "500M"));
        let a = assemble(&["assemble", "-r", "r.fq", "--target-coverage", "20", "--minimap-batch-size", "100M"]);
        assert_eq!((a.reads.target_coverage, a.reads.minimap_batch_size.as_str()), (20, "100M"));
    }

    #[test]
    fn rescue_is_enabled_by_default_and_can_be_disabled() {
        assert!(assemble(&["assemble", "-r", "r.fq"]).rescue_enabled());
        assert!(!assemble(&["assemble", "-r", "r.fq", "--no-rescue"]).rescue_enabled());
        assert!(parse(&["assemble", "-r", "r.fq", "--rescue"]).is_err());
    }
}