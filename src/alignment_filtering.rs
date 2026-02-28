/// alignment filtering module
/// runs in three phases:
/// 1) Read all alignments, store them if they pass basic filters (no self-alignment, overlap length, identity), only store the longest alignment per read pair
/// 2) Calculate coverage statistics per reads
/// 3) Classify alignments into internal matches, contained reads, proper overlaps
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};

// enable serialization for debugging purposes
use serde::{Deserialize, Serialize};
use std::io::BufWriter;

/// Struct to hold a read
struct Read {
    id: usize,
    name: String,
    length: u32,
    per_base_coverage: Vec<u32>,
    coverage_start: u32,
    coverage_end: u32,
}

/// Struct to hold an alignment
#[derive(Clone)]
struct Alignment {
    query_name: String,
    query_length: u32,
    query_start: i64,
    query_end: i64,
    strand: char,
    target_name: String,
    target_length: u32,
    target_start: i64,
    target_end: i64,
    num_matching: u32,
    alignment_block_length: u32,
    mapq: u8,
}

// Implement methods for Alignment (kind of like class methods)
impl Alignment {
    fn from_line(line: &str) -> Option<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            return None;
        }

        Some(Self {
            query_name: fields[0].to_string(),
            query_length: fields[1].parse().ok()?,
            query_start: fields[2].parse::<i64>().ok()?,
            query_end: fields[3].parse::<i64>().ok()?,
            strand: fields[4].chars().next().unwrap_or('+'),
            target_name: fields[5].to_string(),
            target_length: fields[6].parse::<u32>().ok()?,
            target_start: fields[7].parse::<i64>().ok()?,
            target_end: fields[8].parse::<i64>().ok()?,
            num_matching: fields[9].parse().ok()?,
            alignment_block_length: fields[10].parse().ok()?,
            mapq: fields[11].parse().ok()?,
        })
    }

    fn is_self_alignment(&self) -> bool {
        self.query_name == self.target_name
    }

    fn percent_identity(&self) -> f32 {
        (self.num_matching as f32 / self.alignment_block_length as f32) * 100.0
    }
}

/// Enum for alignment classification
enum AlignmentType {
    Filtered,
    InternalMatch,
    FirstContained,
    SecondContained,
    ProperOverlap,
}

// Struct to hold overlap information
#[derive(Serialize, Deserialize)]
pub struct Overlap {
    pub source_name: String,
    pub rc_sink_name: String,
    pub sink_name: String,
    pub rc_source_name: String,
    pub edge_len: u32,
    pub rc_edge_len: u32,
    pub overlap_len: u32,
    pub identity: f64,
}

pub struct AlignmentFilteringOutput {
    pub overlaps: HashMap<(usize, usize), Overlap>,
}

impl AlignmentFilteringOutput {
    // serialize the overlaps
    pub fn serialize_overlaps(&self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let file = File::create(path)?;
        let writer = BufWriter::new(file);
        bincode::serialize_into(writer, &self.overlaps)?;

        Ok(())
    }
}

/// Classify alignment and update contained reads set
fn classify_alignment(
    r: &Alignment,
    query_id: usize,
    target_id: usize,
    overlaps: &mut HashMap<(usize, usize), Overlap>,
    _max_overhang: u32,
    overhang_ratio: f64,
    reads: &Vec<Read>,
    min_overlap_length: u32,
) -> AlignmentType {
    // overlaps are a subset of alignments where (in theory) two read edges, one from each read, are part of the alignment
    // this function tries to differentiate between proper overlaps, internal matches, and containments

    // change overlap coordinates based on read coverage
    let query_start = std::cmp::max(r.query_start, reads[query_id].coverage_start as i64);
    let query_end = std::cmp::min(r.query_end, reads[query_id].coverage_end as i64);
    let target_start = std::cmp::max(r.target_start, reads[target_id].coverage_start as i64);
    let target_end = std::cmp::min(r.target_end, reads[target_id].coverage_end as i64);
    let query_length = (reads[query_id].coverage_end - reads[query_id].coverage_start) as i64;
    let target_length = (reads[target_id].coverage_end - reads[target_id].coverage_start) as i64;

    // using naming convention corresponding with miniasm paper
    let b1 = query_start;
    let e1 = query_end;
    let l1 = query_length;

    // define overlap beginning and end based on orientation
    let (b2, e2, l2) = if r.strand == '+' {
        (target_start, target_end, target_length)
    } else {
        // reverse complement coordinates on the target
        (
            target_length - target_end,
            target_length - target_start,
            target_length,
        )
    };

    // overhang is the part next to the overlap where the reads don't align, but should in case of perfect overlap
    // overhang: min(b1,b2) + min(l1 - e1, l2 - e2)
    let overhang_left = std::cmp::min(b1, b2);
    let overhang_right = std::cmp::min(l1 - e1, l2 - e2);
    let overhang = overhang_left + overhang_right;

    // miniasm uses a hard cutoff on max overhang, we don't for now
    //if overhang_left > max_overhang as i64 || overhang_right > max_overhang as i64 {
    //    return AlignmentType::InternalMatch;
    //}

    // longest overlap length (max aligned spans on either read)
    let overlap_length1 = e1 - b1;
    let overlap_length2 = e2 - b2;
    let overlap_length = std::cmp::max(overlap_length1, overlap_length2) as f64;

    // decide overhang threshold: max_overhang or maplen * overhang_ratio
    // we only use the ratio for now
    let overhang_threshold = (overlap_length * overhang_ratio).ceil() as i64;
    //let max_overhang_i64 = max_overhang as i64;
    //let allowed_overhang = std::cmp::min(max_overhang_i64, overhang_threshold);
    let allowed_overhang = overhang_threshold;

    // classification of the overlap

    // filter out internal matches
    if overhang > allowed_overhang {
        // internal match
        return AlignmentType::InternalMatch;
    }

    // conditions for containment:
    // first contained in second:
    let first_contained = (b1 <= b2) && ((l1 - e1) <= (l2 - e2));
    // second contained in first:
    let second_contained = (b1 >= b2) && ((l1 - e1) >= (l2 - e2));

    // filter out containments
    if first_contained {
        return AlignmentType::FirstContained;
    } else if second_contained {
        return AlignmentType::SecondContained;
    }

    // filter out alignments with very small overlaps
    if overlap_length1 + overhang_left + overhang_right < min_overlap_length as i64
        || overlap_length2 + overhang_left + overhang_right < min_overlap_length as i64
    {
        return AlignmentType::Filtered;
    }

    // at this point it's a proper overlap between reads
    // decide orientation and edge lengths
    if b1 > b2 {
        // first read to second read overlap (query -> target)

        // get nodes & edges (forward orientation)
        let q_plus = format!(
            "{}:{}-{}+",
            r.query_name, reads[query_id].coverage_start, reads[query_id].coverage_end
        );
        let t_orient = format!(
            "{}:{}-{}{}",
            r.target_name, reads[target_id].coverage_start, reads[target_id].coverage_end, r.strand
        );

        // edge length = b1 - b2 (non-overlapping prefix length)
        let edge1_len_i64 = b1 - b2;
        let edge1_len = edge1_len_i64 as u32;

        // reverse complement counterpart:
        // direction: t_rc -> q_minus
        let q_minus = format!(
            "{}:{}-{}-",
            r.query_name, reads[query_id].coverage_start, reads[query_id].coverage_end
        );
        let rc_strand = if r.strand == '+' { '-' } else { '+' };
        let t_rc = format!(
            "{}:{}-{}{}",
            r.target_name,
            reads[target_id].coverage_start,
            reads[target_id].coverage_end,
            rc_strand
        );

        // edge length = (l2 - e2) - (l1 - e1)
        let edge2_len_i64 = (l2 - e2) - (l1 - e1);
        let edge2_len = edge2_len_i64 as u32;

        // shared stats
        let _overlap_len = r.alignment_block_length;
        let identity = r.percent_identity() as f64;

        // create overlap object and store if the overlap is valid
        let ov = Overlap {
            source_name: q_plus,
            rc_sink_name: q_minus,
            sink_name: t_orient,
            rc_source_name: t_rc,
            edge_len: edge1_len,
            rc_edge_len: edge2_len,
            overlap_len: overlap_length as u32,
            identity,
        };
        overlaps.insert((query_id, target_id), ov);
    } else {
        // second to first overlap (target -> query)
        // direction t -> q
        let q_plus = format!(
            "{}:{}-{}+",
            r.query_name, reads[query_id].coverage_start, reads[query_id].coverage_end
        );
        let t_orient = format!(
            "{}:{}-{}{}",
            r.target_name, reads[target_id].coverage_start, reads[target_id].coverage_end, r.strand
        );

        let edge1_len_i64 = b2 - b1;
        let edge1_len = edge1_len_i64 as u32;

        // reverse complement counterpart:
        // direction q_minus -> t_rc
        let q_minus = format!(
            "{}:{}-{}-",
            r.query_name, reads[query_id].coverage_start, reads[query_id].coverage_end
        );
        let rc_strand = if r.strand == '+' { '-' } else { '+' };
        let t_rc = format!(
            "{}:{}-{}{}",
            r.target_name,
            reads[target_id].coverage_start,
            reads[target_id].coverage_end,
            rc_strand
        );

        let edge2_len_i64 = (l1 - e1) - (l2 - e2);
        let edge2_len = edge2_len_i64 as u32;

        // shared stats
        let _overlap_len = r.alignment_block_length;
        let identity = r.percent_identity() as f64;

        // create overlap object and store
        let ov = Overlap {
            source_name: t_orient,
            rc_source_name: q_minus,
            sink_name: q_plus,
            rc_sink_name: t_rc,
            edge_len: edge1_len,
            rc_edge_len: edge2_len,
            overlap_len: overlap_length as u32,
            identity,
        };
        overlaps.insert((query_id, target_id), ov);
    }
    AlignmentType::ProperOverlap
}

/// Filter PAF file based on overlap quality criteria
pub fn run_alignment_filtering(
    paf_in: &str,
    min_overlap_length: &u32,
    min_overlap_count: &u32,
    min_percent_identity: &f32,
    overhang_ratio: &f32,
) -> Result<AlignmentFilteringOutput, Box<dyn std::error::Error>> {
    // Setup data structures
    // read name to read id mapping
    let mut read_name2read_id: HashMap<String, usize> = HashMap::new();
    // vector of read objects, indexed by read id
    let mut reads: Vec<Read> = Vec::new();
    // maps tuples of read ids to alignment objects
    let mut alignments: HashMap<(usize, usize), Alignment> = HashMap::new();
    // vector to track existing alignments per read id
    // useful for querying alignments
    let mut alignment_ids_per_read: Vec<HashSet<usize>> = Vec::new();
    // initialize read id
    let mut next_id: usize = 0;
    // keep track of contained reads
    let mut contained_reads: HashSet<usize> = HashSet::new();
    // initialize overlap storage
    let mut overlaps: HashMap<(usize, usize), Overlap> = HashMap::new();

    let mut self_alignments_skipped: usize = 0;
    let mut alignment_length_skipped: usize = 0;
    let mut percent_identity_skipped: usize = 0;

    // read the alignments from the PAF file
    let reader = BufReader::new(File::open(paf_in)?);
    for line in reader.lines() {
        let line = line?;

        // skip header lines
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        if let Some(record) = Alignment::from_line(&line) {
            // skip self alignments
            if record.is_self_alignment() {
                self_alignments_skipped += 1;
                continue;
            }

            let query_overlap_length = record.query_end - record.query_start;
            let target_overlap_length = record.target_end - record.target_start;

            // skip low overlap length alignments
            if query_overlap_length < (*min_overlap_length).into()
                || target_overlap_length < (*min_overlap_length).into()
            {
                alignment_length_skipped += 1;
                continue;
            }

            // skip low percent identity alignments
            if record.percent_identity() < *min_percent_identity {
                percent_identity_skipped += 1;
                continue;
            }
            // Get or create read ids for query and target
            let query_id = match read_name2read_id.get(&record.query_name) {
                // we have seen this read before, get its id
                Some(&id) => id,

                // new read, assign new id and create read object
                None => {
                    let id = next_id;
                    next_id += 1;
                    read_name2read_id.insert(record.query_name.clone(), id);
                    alignment_ids_per_read.push(HashSet::new());
                    // create new read object
                    reads.push(Read {
                        id,
                        name: record.query_name.clone(),
                        length: record.query_length,
                        per_base_coverage: vec![0; record.query_length as usize],
                        coverage_start: 0,
                        coverage_end: record.query_length,
                    });
                    id
                }
            };

            let target_id = match read_name2read_id.get(&record.target_name) {
                Some(&id) => id,
                None => {
                    let id = next_id;
                    next_id += 1;
                    read_name2read_id.insert(record.target_name.clone(), id);
                    alignment_ids_per_read.push(HashSet::new());
                    // create new read object
                    reads.push(Read {
                        id,
                        name: record.target_name.clone(),
                        length: record.target_length,
                        per_base_coverage: vec![0; record.target_length as usize],
                        coverage_start: 0,
                        coverage_end: record.target_length,
                    });
                    id
                }
            };

            // extract needed info from the record
            let (qstart, qend) = (record.query_start as usize, record.query_end as usize);
            let (tstart, tend) = (record.target_start as usize, record.target_end as usize);

            // store alignment record
            // if multiple alignments exist between the same read pair, keep the longest one

            if alignment_ids_per_read[query_id].contains(&target_id) {
                // an alignment between these reads already exists, it may be stored under
                // (query_id, target_id) or (target_id, query_id)
                if let Some(existing) = alignments.get(&(query_id, target_id)) {
                    if record.alignment_block_length > existing.alignment_block_length {
                        // replace the existing directed entry
                        alignments.insert((query_id, target_id), record.clone());
                    }
                } else if let Some(existing) = alignments.get(&(target_id, query_id)) {
                    if record.alignment_block_length > existing.alignment_block_length {
                        // remove the old reversed entry and store the new (keeps orientation of current record)
                        alignments.remove(&(target_id, query_id));
                        alignments.insert((query_id, target_id), record.clone());
                    }
                } else {
                    println!("Warning: alignment existence inconsistency detected.");
                }
            }
            // we don't have an alignment between these reads yet
            else {
                alignment_ids_per_read[query_id].insert(target_id);
                alignment_ids_per_read[target_id].insert(query_id);
                alignments.insert((query_id, target_id), record);
            }

            // update read coverage statistics
            reads[query_id].per_base_coverage[qstart..qend]
                .iter_mut()
                .for_each(|c| *c += 1);
            reads[target_id].per_base_coverage[tstart..tend]
                .iter_mut()
                .for_each(|c| *c += 1);
        }
    }

    println!("=== ALIGNMENT FILTERING ===");
    println!("=== PHASE 1: CRUDE FILTERING ===");
    println!("Total self-alignments skipped: {}", self_alignments_skipped);
    println!(
        "Total alignments skipped due to length filter: {}",
        alignment_length_skipped
    );
    println!(
        "Total alignments skipped due to percent identity filter: {}",
        percent_identity_skipped
    );
    println!("Total reads kept: {}", reads.len());
    println!("Total alignments kept: {}", alignments.len());
    println!("=== PHASE 1 FINISHED ===");
    println!("=== PHASE 2: COVERAGE CALCULATION ===");

    // all alignments have been read
    // store subregions with coverage >= 3 per read
    for read in &mut reads {
        let mut cur_len = 0;
        let mut cur_start = 0;
        let mut best_len = 0;
        let mut best = (0, 0);
        for (i, x) in read.per_base_coverage.iter().enumerate() {
            if x > min_overlap_count {
                if cur_len == 0 {
                    cur_start = i;
                }
                cur_len += 1;
                if cur_len > best_len {
                    best_len = cur_len;
                    best = (cur_start, i + 1); // end is exclusive
                }
            } else {
                cur_len = 0;
            }
        }
        read.coverage_start = best.0 as u32;
        read.coverage_end = best.1 as u32;
    }

    println!("=== PHASE 2 FINISHED ===");
    println!("=== PHASE 3: ALIGNMENT CLASSIFICATION ===");

    // classify alignments and update contained reads set
    for ((query_id, target_id), alignment) in &alignments {
        //let alignment_type = match classify_alignment(alignment, *query_id, *target_id, &mut overlaps, *max_overhang, *overhang_ratio, &reads, min_overlap_length) {
        let _alignment_type = match classify_alignment(
            alignment,
            *query_id,
            *target_id,
            &mut overlaps,
            4000000000,
            (*overhang_ratio) as f64,
            &reads,
            *min_overlap_length,
        ) {
            AlignmentType::Filtered => {
                continue; // skip filtered alignments
            }
            AlignmentType::InternalMatch => {
                continue; // skip internal matches
            }
            AlignmentType::FirstContained => {
                // mark first read as contained
                contained_reads.insert(*query_id);
                continue;
            }
            AlignmentType::SecondContained => {
                // mark second read as contained
                contained_reads.insert(*target_id);
                continue;
            }
            AlignmentType::ProperOverlap => {
                continue;
            }
        };
    }

    println!("Total overlaps after classification: {}", overlaps.len());

    // filter contained reads from overlaps
    overlaps.retain(|(q_id, t_id), _| {
        !contained_reads.contains(q_id) && !contained_reads.contains(t_id)
    });

    println!(
        "Total overlaps after removing contained reads: {}",
        overlaps.len()
    );

    // get unique reads from overlaps
    let unique_reads: HashSet<usize> = overlaps
        .keys()
        .flat_map(|(q_id, t_id)| vec![*q_id, *t_id])
        .collect();

    // filter low coverage reads
    let threshold = *min_overlap_count;
    let low_coverage_reads: Vec<_> = reads
        .iter()
        .enumerate()
        .filter(|(_, r)| r.per_base_coverage.iter().copied().max().unwrap_or(0) < threshold)
        .map(|(id, _)| id)
        .collect();
    overlaps.retain(|(q_id, t_id), _| {
        !low_coverage_reads.contains(q_id) && !low_coverage_reads.contains(t_id)
    });
    println!(
        "Total number of reads for graph creation: {}",
        unique_reads.len()
    );
    println!(
        "Total number of overlaps for graph creation: {}",
        overlaps.len()
    );
    println!("=== PHASE 3 FINISHED ===");
    println!("=== ALIGNMENT FILTERING FINISHED ===");

    Ok(AlignmentFilteringOutput { overlaps })
}
