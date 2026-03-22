pub struct AlignReadsConfig {
    // output
    pub output_dir: String,

    // read filtering and alignment
    pub reads_fq: String,
    pub threads: usize,
    pub paf: String,
    pub min_read_length: u32,
    pub min_base_quality: f32,
}

pub struct AlignmentFilteringConfig {
    pub output_dir: String,
    pub paf: String,
    pub output_overlaps: String,
    pub min_overlap_length: u32,
    pub min_overlap_count: u32,
    pub min_percent_identity: f32,
    pub overhang_ratio: f32,
}

pub struct AssembleConfig {
    // output
    pub output_prefix: String,
    pub output_dir: String,

    // read filtering and alignment
    pub reads_fq: String,
    pub threads: usize,
    pub paf: String,
    pub min_read_length: u32,
    pub min_base_quality: f32,

    // alignment filtering
    pub min_overlap_length: u32,
    pub min_overlap_count: u32,
    pub min_percent_identity: f32,
    pub overhang_ratio: f32,
    pub overlaps: Option<String>,
    
    // assembly parameters
    pub max_bubble_length: u32,
    pub min_support_ratio: f64,
    pub max_tip_len: u32,
    pub fuzz: u32,
    pub cleanup_iterations: u32,
    pub short_edge_ratio: f64,
}
