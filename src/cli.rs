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
    version = "1.2.1",
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
}

#[derive(Args)]
pub struct AlignReadsArgs {
    /// Output directory
    #[arg(short = 'o', long, default_value = ".")]
    pub output_dir: String,

    /// Input reads in FASTQ format
    #[arg(short = 'r', long)]
    pub reads_fq: String,

    /// Read type / sequencing technology
    #[arg(long, default_value = "ont")]
    pub read_type: ReadType,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = 4)]
    pub threads: usize,

    /// Output PAF file
    #[arg(short = 'a', long, default_value = "alignments.paf")]
    pub paf: String,

    /// Minimum read length
    #[arg(long, default_value_t = 1000)]
    pub min_read_length: u32,

    /// Minimum average quality
    #[arg(short = 'q', long, default_value_t = 10.0)]
    pub min_base_quality: f32,

    /// Optional input genome size (if not provided, will be estimated from data)
    #[arg(long)]
    pub genome_size: Option<u32>,

    /// Target read coverage for assembly input downsampling
    #[arg(long, default_value_t = 50u32)]
    pub target_coverage: u32,

    /// Minimap2 query minibatch size for self-alignment (e.g. 500M)
    #[arg(long, default_value = "500M")]
    pub minimap_batch_size: String,

    /// Override minimap2's -f high-frequency minimizer filter
    #[arg(long)]
    pub minimap2_f: Option<String>,

    /// Override minimap2's -U minimizer occurrence bounds
    #[arg(long)]
    pub minimap2_u: Option<String>,
}

impl From<&AlignReadsArgs> for crate::configs::AlignReadsConfig {
    fn from(args: &AlignReadsArgs) -> Self {
        Self {
            output_dir: args.output_dir.clone(),
            reads_fq: args.reads_fq.clone(),
            threads: args.threads,
            paf: args.paf.clone(),
            min_read_length: args.min_read_length,
            min_base_quality: args.min_base_quality,
            genome_size: args.genome_size,
            target_coverage: args.target_coverage,
            minimap_batch_size: args.minimap_batch_size.clone(),
            minimap2_f: args.minimap2_f.clone(),
            minimap2_u: args.minimap2_u.clone(),
            read_type: args.read_type,
        }
    }
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

    /// Optional random seed for reproducible results
    #[arg(long)]
    pub seed: Option<u64>,
}

impl From<&AlignmentFilteringArgs> for crate::configs::AlignmentFilteringConfig {
    fn from(args: &AlignmentFilteringArgs) -> Self {
        Self {
            output_dir: args.output_dir.clone(),
            paf: args.paf.clone(),
            output_overlaps: args.output_overlaps.clone(),
            min_overlap_length: args.min_overlap_length,
            min_overlap_count: args.min_overlap_count,
            min_percent_identity: args.min_percent_identity,
            overhang_ratio: args.overhang_ratio,
            seed: args.seed,
        }
    }
}

#[derive(Args)]
pub struct AssembleArgs {
    /// Output parameters

    /// Output prefix
    #[arg(short = 'p', long, default_value = "unitigs", help_heading = "Output")]
    pub output_prefix: String,

    /// Output directory
    #[arg(short = 'o', long, default_value = ".", help_heading = "Output")]
    pub output_dir: String,

    /// Read filtering and alignment parameters

    /// Input reads in FASTQ format
    #[arg(short = 'r', long, help_heading = "Read filtering and alignment")]
    pub reads_fq: String,

    /// Read type / sequencing technology
    #[arg(
        long,
        default_value = "ont",
        help_heading = "Read filtering and alignment"
    )]
    pub read_type: ReadType,

    /// Number of threads
    #[arg(
        short = 't',
        long,
        default_value_t = 4,
        help_heading = "Read filtering and alignment"
    )]
    pub threads: usize,

    /// Output PAF filename
    #[arg(
        short = 'a',
        long,
        default_value = "alignments.paf",
        help_heading = "Read filtering and alignment"
    )]
    pub paf: String,

    /// Minimum read length
    #[arg(
        long,
        default_value_t = 1000,
        help_heading = "Read filtering and alignment"
    )]
    pub min_read_length: u32,

    /// Minimum average quality
    #[arg(
        short = 'q',
        long,
        default_value_t = 10.0,
        help_heading = "Read filtering and alignment"
    )]
    pub min_base_quality: f32,

    /// Optional input genome size (if not provided, will be estimated from data)
    #[arg(long, help_heading = "Read filtering and alignment")]
    pub genome_size: Option<u32>,

    /// Target read coverage for assembly input downsampling
    #[arg(
        long,
        default_value_t = 50u32,
        help_heading = "Read filtering and alignment"
    )]
    pub target_coverage: u32,

    /// Minimap2 query minibatch size for self-alignment (e.g. 500M)
    #[arg(
        long,
        default_value = "500M",
        help_heading = "Read filtering and alignment"
    )]
    pub minimap_batch_size: String,

    /// Override minimap2's -f high-frequency minimizer filter
    #[arg(long, help_heading = "Read filtering and alignment")]
    pub minimap2_f: Option<String>,

    /// Override minimap2's -U minimizer occurrence bounds
    #[arg(long, help_heading = "Read filtering and alignment")]
    pub minimap2_u: Option<String>,

    /// Alignment filtering parameters (optional if --overlaps is provided)

    /// Minimum overlap length
    #[arg(
        short = 'l',
        long,
        default_value_t = 2000,
        help_heading = "Alignment filtering"
    )]
    pub min_overlap_length: u32,

    /// Minimum overlap count
    #[arg(
        short = 'c',
        long,
        default_value_t = 3,
        help_heading = "Alignment filtering"
    )]
    pub min_overlap_count: u32,

    /// Minimum percent identity
    #[arg(
        short = 'i',
        long,
        default_value_t = 5.0,
        help_heading = "Alignment filtering"
    )]
    pub min_percent_identity: f32,

    /// Overhang ratio
    #[arg(long, default_value_t = 0.8, help_heading = "Alignment filtering")]
    pub overhang_ratio: f32,

    /// Pre-computed overlaps binary file (optional, if provided skips alignment filtering)
    #[arg(long, help_heading = "Alignment filtering")]
    pub overlaps: Option<String>,

    /// Assembly parameters

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

    /// Enable a second read-guided completion stage after the first assembly pass
    #[arg(long, default_value_t = false, action = clap::ArgAction::SetTrue, help_heading = "Assembly")]
    pub completion_enabled: bool,

    /// Number of completion rounds to run
    #[arg(long, default_value_t = 1u32, help_heading = "Assembly")]
    pub completion_rounds: u32,

    /// Minimum alignment length for a completion bridge
    #[arg(long, default_value_t = 2000u32, help_heading = "Assembly")]
    pub completion_min_alignment_len: u32,

    /// Minimum identity for a completion bridge
    #[arg(long, default_value_t = 0.8f64, help_heading = "Assembly")]
    pub completion_min_identity: f64,

    /// Assemble poorly represented reads once as a plasmid rescue stage
    #[arg(long, default_value_t = true, action = clap::ArgAction::SetTrue, help_heading = "Assembly")]
    pub rescue_plasmids: bool,

    #[arg(long, hide = true, action = clap::ArgAction::SetTrue)]
    pub no_rescue_plasmids: bool,

    /// Optional random seed for reproducible results
    #[arg(long, help_heading = "Assembly")]
    pub seed: Option<u64>,
}

impl From<&AssembleArgs> for crate::configs::AssembleConfig {
    fn from(args: &AssembleArgs) -> Self {
        Self {
            // output
            output_prefix: args.output_prefix.clone(),
            output_dir: args.output_dir.clone(),

            // read filtering and alignment
            reads_fq: args.reads_fq.clone(),
            threads: args.threads,
            paf: args.paf.clone(),
            min_read_length: args.min_read_length,
            min_base_quality: args.min_base_quality,
            genome_size: args.genome_size,
            target_coverage: args.target_coverage,
            minimap_batch_size: args.minimap_batch_size.clone(),
            minimap2_f: args.minimap2_f.clone(),
            minimap2_u: args.minimap2_u.clone(),
            read_type: args.read_type,

            // alignment filtering
            min_overlap_length: args.min_overlap_length,
            min_overlap_count: args.min_overlap_count,
            min_percent_identity: args.min_percent_identity,
            overhang_ratio: args.overhang_ratio,
            overlaps: args.overlaps.clone(),

            // assembly parameters
            max_bubble_length: args.max_bubble_length,
            min_support_ratio: args.min_support_ratio,
            max_tip_len: args.max_tip_len,
            fuzz: args.fuzz,
            cleanup_iterations: args.cleanup_iterations,
            short_edge_ratio: args.short_edge_ratio,
            completion_enabled: args.completion_enabled,
            completion_rounds: args.completion_rounds,
            completion_min_alignment_len: args.completion_min_alignment_len,
            completion_min_identity: args.completion_min_identity,
            rescue_plasmids: args.rescue_plasmids && !args.no_rescue_plasmids,
            seed: args.seed,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_minimap2_args_for_read_types() {
        assert_eq!(ReadType::Ont.minimap2_args(), &["-x", "ava-ont"]);
        assert_eq!(ReadType::PbClr.minimap2_args(), &["-x", "ava-pb"]);
        assert_eq!(
            ReadType::PbHifi.minimap2_args(),
            &[
                "-x", "ava-ont", "-k", "21", "-w", "11", "-g", "1000", "-m", "200", "-r", "2000",
            ]
        );
        assert_eq!(ReadType::Ont.mapping_args(), &["-x", "map-ont"]);
        assert_eq!(ReadType::PbClr.mapping_args(), &["-x", "map-pb"]);
        assert_eq!(ReadType::PbHifi.mapping_args(), &["-x", "map-hifi"]);
    }

    #[test]
    fn test_cli_read_type_default() {
        let cli = Cli::try_parse_from(["Ilesta", "align", "-r", "reads.fq"]).unwrap();
        match cli.command {
            Commands::Align(args) => assert_eq!(args.read_type, ReadType::Ont),
            _ => panic!("Expected Align command"),
        }

        let cli_asm = Cli::try_parse_from(["Ilesta", "assemble", "-r", "reads.fq"]).unwrap();
        match cli_asm.command {
            Commands::Assemble(args) => assert_eq!(args.read_type, ReadType::Ont),
            _ => panic!("Expected Assemble command"),
        }
    }

    #[test]
    fn test_cli_read_type_explicit() {
        for (flag_val, expected) in [
            ("ont", ReadType::Ont),
            ("pb-clr", ReadType::PbClr),
            ("pb-hifi", ReadType::PbHifi),
        ] {
            let cli =
                Cli::try_parse_from(["Ilesta", "align", "-r", "reads.fq", "--read-type", flag_val])
                    .unwrap();
            match cli.command {
                Commands::Align(args) => assert_eq!(args.read_type, expected),
                _ => panic!("Expected Align command"),
            }

            let cli_asm = Cli::try_parse_from([
                "Ilesta",
                "assemble",
                "-r",
                "reads.fq",
                "--read-type",
                flag_val,
            ])
            .unwrap();
            match cli_asm.command {
                Commands::Assemble(args) => assert_eq!(args.read_type, expected),
                _ => panic!("Expected Assemble command"),
            }
        }
    }

    #[test]
    fn test_cli_read_type_invalid() {
        let result = Cli::try_parse_from([
            "Ilesta",
            "align",
            "-r",
            "reads.fq",
            "--read-type",
            "illumina",
        ]);
        assert!(result.is_err());
    }

    #[test]
    fn test_align_target_coverage_option() {
        let default_cli = Cli::try_parse_from(["Ilesta", "align", "-r", "reads.fq"]).unwrap();
        match default_cli.command {
            Commands::Align(args) => assert_eq!(args.target_coverage, 50),
            _ => panic!("Expected Align command"),
        }

        let configured_cli = Cli::try_parse_from([
            "Ilesta",
            "align",
            "-r",
            "reads.fq",
            "--target-coverage",
            "20",
        ])
        .unwrap();
        match configured_cli.command {
            Commands::Align(args) => assert_eq!(args.target_coverage, 20),
            _ => panic!("Expected Align command"),
        }
    }

    #[test]
    fn test_minimap_batch_size_defaults_and_overrides() {
        let default_align = Cli::try_parse_from(["Ilesta", "align", "-r", "reads.fq"]).unwrap();
        match default_align.command {
            Commands::Align(args) => assert_eq!(args.minimap_batch_size, "500M"),
            _ => panic!("Expected Align command"),
        }

        let configured_assemble = Cli::try_parse_from([
            "Ilesta",
            "assemble",
            "-r",
            "reads.fq",
            "--minimap-batch-size",
            "100M",
        ])
        .unwrap();
        match configured_assemble.command {
            Commands::Assemble(args) => assert_eq!(args.minimap_batch_size, "100M"),
            _ => panic!("Expected Assemble command"),
        }
    }

    #[test]
    fn test_cli_rescue_plasmids_flag() {
        let cli =
            Cli::try_parse_from(["Ilesta", "assemble", "-r", "reads.fq", "--rescue-plasmids"])
                .unwrap();
        match cli.command {
            Commands::Assemble(args) => assert!(args.rescue_plasmids),
            _ => panic!("Expected Assemble command"),
        }
    }

    #[test]
    fn test_cli_can_disable_nested_plasmid_rescue() {
        let cli = Cli::try_parse_from([
            "Ilesta",
            "assemble",
            "-r",
            "reads.fq",
            "--no-rescue-plasmids",
        ])
        .unwrap();
        match cli.command {
            Commands::Assemble(args) => {
                let config: crate::configs::AssembleConfig = (&args).into();
                assert!(!config.rescue_plasmids);
            }
            _ => panic!("Expected Assemble command"),
        }
    }
}
