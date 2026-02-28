use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[command(
    name = "Ilesta",
    version = "1.1.0",
    about = "De novo genome assembly for long reads using an overlap graph"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Alignment filtering
    AlignmentFiltering(AlignmentFilteringArgs),

    /// Full genome assembly pipeline
    Assemble(AssembleArgs),
}

#[derive(Args)]
pub struct AlignmentFilteringArgs {
    /// Input PAF file
    #[arg(short = 'f', long)]
    pub input_paf: String,

    /// Output overlaps binary file
    #[arg(long, default_value = "overlaps.bin")]
    pub output_overlaps: String,

    /// Minimum overlap length
    #[arg(short = 'l', long, default_value_t = 2000)]
    pub min_overlap_length: u32,

    /// Minimum overlap count
    #[arg(short = 'c', long, default_value_t = 3)]
    pub min_overlap_count: u32,

    /// Minimum percent identity
    #[arg(short = 'i', long, default_value_t = 5.0)]
    pub min_percent_identity: f32,

    /// Overhang ratio
    #[arg(long, default_value_t = 0.8)]
    pub overhang_ratio: f32,
}

impl From<&AlignmentFilteringArgs> for crate::configs::AlignmentFilteringConfig {
    fn from(args: &AlignmentFilteringArgs) -> Self {
        Self {
            input_paf: args.input_paf.clone(),
            output_overlaps: args.output_overlaps.clone(),
            min_overlap_length: args.min_overlap_length,
            min_overlap_count: args.min_overlap_count,
            min_percent_identity: args.min_percent_identity,
            overhang_ratio: args.overhang_ratio,
        }
    }
}

#[derive(Args)]
pub struct AssembleArgs {
    // Alignment filtering parameters (optional if --overlaps is provided)
    /// Input PAF file (optional if --overlaps is provided)
    #[arg(short = 'f', long)]
    pub input_paf: Option<String>,

    /// Minimum overlap length
    #[arg(short = 'l', long, default_value_t = 2000)]
    pub min_overlap_length: u32,

    /// Minimum overlap count
    #[arg(short = 'c', long, default_value_t = 3)]
    pub min_overlap_count: u32,

    /// Minimum percent identity
    #[arg(short = 'i', long, default_value_t = 5.0)]
    pub min_percent_identity: f32,

    /// Overhang ratio
    #[arg(long, default_value_t = 0.8)]
    pub overhang_ratio: f32,

    /// Pre-computed overlaps binary file (optional, if provided skips alignment filtering)
    #[arg(long)]
    pub overlaps: Option<String>,

    /// Input reads in FASTQ format
    #[arg(short = 'r', long)]
    pub reads_fq: String,

    /// Output prefix
    #[arg(short = 'p', long, default_value = "unitigs")]
    pub output_prefix: String,

    /// Output directory
    #[arg(short = 'o', long, default_value = ".")]
    pub output_dir: String,

    /// Maximum bubble length (used during bubble removal)
    #[arg(long, default_value_t = 100u32)]
    pub max_bubble_length: u32,

    /// Minimum support ratio for bubble removal
    #[arg(long, default_value_t = 1.1f64)]
    pub min_support_ratio: f64,

    /// Maximum tip length for tip trimming
    #[arg(long, default_value_t = 4u32)]
    pub max_tip_len: u32,

    /// Fuzz parameter for transitive edge reduction
    #[arg(long, default_value_t = 10u32)]
    pub fuzz: u32,

    /// Number of cleanup iterations to run
    #[arg(long, default_value_t = 3u32)]
    pub cleanup_iterations: u32,

    /// Short edge removal ratio (heuristic simplification)
    #[arg(long, default_value_t = 0.8f64)]
    pub short_edge_ratio: f64,
}

impl From<&AssembleArgs> for crate::configs::AssembleConfig {
    fn from(args: &AssembleArgs) -> Self {
        Self {
            input_paf: args.input_paf.clone(),
            min_overlap_length: args.min_overlap_length,
            min_overlap_count: args.min_overlap_count,
            min_percent_identity: args.min_percent_identity,
            overhang_ratio: args.overhang_ratio,
            overlaps: args.overlaps.clone(),
            reads_fq: args.reads_fq.clone(),
            output_prefix: args.output_prefix.clone(),
            output_dir: args.output_dir.clone(),
            max_bubble_length: args.max_bubble_length,
            min_support_ratio: args.min_support_ratio,
            max_tip_len: args.max_tip_len,
            fuzz: args.fuzz,
            cleanup_iterations: args.cleanup_iterations,
            short_edge_ratio: args.short_edge_ratio,
        }
    }
}
