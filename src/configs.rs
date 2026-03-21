pub struct AlignReadsConfig {
    pub reads_fq: String,
    pub threads: usize,
    pub output_paf: String,
}

pub struct AlignmentFilteringConfig {
    pub input_paf: String,
    pub output_overlaps: String,
    pub min_overlap_length: u32,
    pub min_overlap_count: u32,
    pub min_percent_identity: f32,
    pub overhang_ratio: f32,
}

pub struct AssembleConfig {
    pub input_paf: Option<String>,
    pub min_overlap_length: u32,
    pub min_overlap_count: u32,
    pub min_percent_identity: f32,
    pub overhang_ratio: f32,
    pub overlaps: Option<String>,
    pub reads_fq: String,
    pub output_prefix: String,
    pub output_dir: String,
    pub max_bubble_length: u32,
    pub min_support_ratio: f64,
    pub max_tip_len: u32,
    pub fuzz: u32,
    pub cleanup_iterations: u32,
    pub short_edge_ratio: f64,
}
