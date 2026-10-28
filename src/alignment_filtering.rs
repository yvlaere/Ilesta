/// alignment filtering module
/// runs in three phases:
/// 1) Read all alignments, store them if they pass basic filters (no self-alignment, overlap length, identity), only store the longest alignment per read pair
/// 2) Calculate coverage statistics per reads
/// 3) Classify alignments into internal matches, contained reads, proper overlaps
use crate::ui;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};

// enable serialization for debugging purposes
use serde::{Deserialize, Serialize};
use std::io::BufWriter;

/// Struct to hold a read
struct Read {
    per_base_coverage: Vec<u16>,
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
        })
    }

    fn is_self_alignment(&self) -> bool {
        self.query_name == self.target_name
    }

    /// A read aligned to its own reverse complement (strand '-') with the two arms adjacent: the
    /// read folds back on itself (a chimeric hairpin read, second half = reverse complement of the
    /// first). Arm length and gap are fractions of the read length, so that short inverted repeats
    /// (IS elements, rRNA operons) in a long read do not qualify. A hit that spans the turning
    /// point is reported as one record (query and target overlap), then the gap is 0.
    fn is_fold_back(&self) -> bool {
        let gap = (self.query_start.max(self.target_start) - self.query_end.min(self.target_end)).max(0);
        let arm = (self.query_end.max(self.target_end) - self.query_start.min(self.target_start) - gap) / 2;
        let len = self.query_length as f64;
        self.is_self_alignment()
            && self.strand == '-'
            && arm as f64 >= MIN_FOLD_BACK_ARM * len
            && gap as f64 <= MAX_FOLD_BACK_GAP * len
    }

    fn percent_identity(&self) -> f32 {
        (self.num_matching as f32 / self.alignment_block_length as f32) * 100.0
    }
}

/// Shortest self-inverted arm of a fold-back read, as a fraction of the read length.
const MIN_FOLD_BACK_ARM: f64 = 0.2;
/// Largest gap between the two arms (hairpin loop, adapter remnants), as a fraction of the read length.
const MAX_FOLD_BACK_GAP: f64 = 0.01;

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
        if crate::utils::has_seed() {
            let sorted: std::collections::BTreeMap<_, _> = self.overlaps.iter().collect();
            bincode::serialize_into(writer, &sorted)?;
        } else {
            bincode::serialize_into(writer, &self.overlaps)?;
        }

        Ok(())
    }
}

pub static CLIPPED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static CLIPPED_BASES: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// Classify alignment and update contained reads set
fn classify_alignment(
    r: &Alignment,
    query_id: usize,
    target_id: usize,
    overlaps: &mut HashMap<(usize, usize), Overlap>,
    overhang_ratio: f64,
    reads: &[Read],
    min_overlap_length: u32,
) -> AlignmentType {
    // overlaps are a subset of alignments where (in theory) two read edges, one from each read, are part of the alignment
    // this function tries to differentiate between proper overlaps, internal matches, and containments

    // change overlap coordinates based on read coverage; all coordinates are relative to the
    // coverage-trimmed read (the node sequence is `read[coverage_start..coverage_end]`)
    let (qcs, qce) = (reads[query_id].coverage_start as i64, reads[query_id].coverage_end as i64);
    let (tcs, tce) = (reads[target_id].coverage_start as i64, reads[target_id].coverage_end as i64);
    // Clipping one read to its trimmed interval must advance the other read by the same number of
    // bases (the alignment diagonal is fixed), otherwise the offset between the two reads, and
    // with it the edge length, is off by the clipped length.
    let plus = r.strand == '+';
    let lead = [qcs - r.query_start, if plus { tcs - r.target_start } else { r.target_end - tce }, 0]
        .into_iter()
        .max()
        .unwrap();
    let trail = [r.query_end - qce, if plus { r.target_end - tce } else { tcs - r.target_start }, 0]
        .into_iter()
        .max()
        .unwrap();
    let (qs, qe) = (r.query_start + lead, r.query_end - trail);
    let (mut ts, mut te) = (r.target_start, r.target_end);
    if plus {
        ts += lead;
        te -= trail;
    } else {
        te -= lead;
        ts += trail;
    }
    if qs >= qe || ts >= te {
        return AlignmentType::Filtered;
    }
    if lead > 0 || trail > 0 {
        CLIPPED.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        CLIPPED_BASES.fetch_add((lead + trail) as u64, std::sync::atomic::Ordering::Relaxed);
    }
    let (query_start, query_end) = (qs - qcs, qe - qcs);
    let (target_start, target_end) = (ts - tcs, te - tcs);
    let query_length = qce - qcs;
    let target_length = tce - tcs;

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

    // longest overlap length (max aligned spans on either read)
    let overlap_length1 = e1 - b1;
    let overlap_length2 = e2 - b2;
    let overlap_length = std::cmp::max(overlap_length1, overlap_length2) as f64;

    // allowed overhang: overlap length * overhang_ratio
    let allowed_overhang = (overlap_length * overhang_ratio).ceil() as i64;

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
    paf_in: &std::path::Path,
    opts: &crate::cli::OverlapOpts,
    report_stats: bool,
) -> Result<AlignmentFilteringOutput, Box<dyn std::error::Error>> {
    let crate::cli::OverlapOpts {
        min_overlap_length,
        min_overlap_count,
        min_percent_identity,
        overhang_ratio,
    } = opts;
    CLIPPED.store(0, std::sync::atomic::Ordering::Relaxed);
    CLIPPED_BASES.store(0, std::sync::atomic::Ordering::Relaxed);
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
    // names of fold-back reads
    let mut fold_back_names: HashMap<String, u32> = HashMap::new();
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
                if record.is_fold_back() {
                    fold_back_names.insert(record.query_name.clone(), record.query_length);
                }
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
                    crate::ui::warn("alignment existence inconsistency detected");
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
                .for_each(|c| *c = c.saturating_add(1));
            reads[target_id].per_base_coverage[tstart..tend]
                .iter_mut()
                .for_each(|c| *c = c.saturating_add(1));
        }
    }

    // all alignments have been read
    // drop every alignment of a fold-back read
    let fold_back_reads: HashSet<usize> =
        fold_back_names.keys().filter_map(|n| read_name2read_id.get(n).copied()).collect();
    if !fold_back_reads.is_empty() {
        alignments.retain(|(q, t), _| !fold_back_reads.contains(q) && !fold_back_reads.contains(t));
    }
    // store subregions with coverage >= 3 per read
    for read in &mut reads {
        let mut cur_len = 0;
        let mut cur_start = 0;
        let mut best_len = 0;
        let mut best = (0, 0);
        for (i, x) in read.per_base_coverage.iter().enumerate() {
            if (*x as u32) > *min_overlap_count {
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

    // classify alignments and update contained reads set
    for ((query_id, target_id), alignment) in &alignments {
        let _alignment_type = match classify_alignment(
            alignment,
            *query_id,
            *target_id,
            &mut overlaps,
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

    // filter contained reads from overlaps
    overlaps.retain(|(q_id, t_id), _| {
        !contained_reads.contains(q_id) && !contained_reads.contains(t_id)
    });

    // filter low coverage reads
    let threshold = *min_overlap_count;
    let low_coverage_reads: HashSet<usize> = reads
        .iter()
        .enumerate()
        .filter(|(_, r)| (r.per_base_coverage.iter().copied().max().unwrap_or(0) as u32) < threshold)
        .map(|(id, _)| id)
        .collect();
    overlaps.retain(|(q_id, t_id), _| {
        !low_coverage_reads.contains(q_id) && !low_coverage_reads.contains(t_id)
    });
    let graph_reads: HashSet<usize> = overlaps
        .keys()
        .flat_map(|(q_id, t_id)| [*q_id, *t_id])
        .collect();
    if report_stats {
        ui::detail(format!(
            "{} reads · {} alignments kept",
            ui::count(reads.len() as u64),
            ui::count(alignments.len() as u64)
        ));
        ui::detail(format!(
            "skipped: {} self, {} short, {} low-identity alignments",
            ui::count(self_alignments_skipped as u64),
            ui::count(alignment_length_skipped as u64),
            ui::count(percent_identity_skipped as u64)
        ));
        ui::detail(format!(
            "fold-back reads removed: {} ({})",
            ui::count(fold_back_names.len() as u64),
            ui::bp(fold_back_names.values().map(|&l| u64::from(l)).sum::<u64>())
        ));
        ui::detail(format!(
            "alignments clipped to coverage-trimmed reads: {} ({})",
            ui::count(CLIPPED.load(std::sync::atomic::Ordering::Relaxed)),
            ui::bp(CLIPPED_BASES.load(std::sync::atomic::Ordering::Relaxed))
        ));
        ui::detail(format!(
            "overlap graph input: {} reads · {} overlaps ({} contained, {} low-coverage reads removed)",
            ui::count(graph_reads.len() as u64),
            ui::count(overlaps.len() as u64),
            ui::count(contained_reads.len() as u64),
            ui::count(low_coverage_reads.len() as u64)
        ));
    }

    Ok(AlignmentFilteringOutput { overlaps })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(start: u32, end: u32) -> Read {
        Read { per_base_coverage: Vec::new(), coverage_start: start, coverage_end: end }
    }

    #[allow(clippy::too_many_arguments)]
    fn aln(qs: i64, qe: i64, strand: char, ts: i64, te: i64) -> Alignment {
        Alignment {
            query_name: "q".into(),
            query_length: 5000,
            query_start: qs,
            query_end: qe,
            strand,
            target_name: "t".into(),
            target_length: 5000,
            target_start: ts,
            target_end: te,
            num_matching: 3900,
            alignment_block_length: 4000,
        }
    }

    fn classify(a: &Alignment, reads: &[Read]) -> Option<Overlap> {
        let mut overlaps = HashMap::new();
        classify_alignment(a, 0, 1, &mut overlaps, 0.2, reads, 500);
        overlaps.remove(&(0, 1))
    }

    /// t starts 1000 bp before q. q is coverage-trimmed to start at 400, so the t -> q edge is
    /// 1000 + 400 bp long (clipping q must also advance t, not leave its start at 1000).
    #[test]
    fn trimmed_query_start_shifts_target() {
        let reads = [read(400, 5000), read(0, 5000)];
        let ov = classify(&aln(0, 4000, '+', 1000, 5000), &reads).expect("overlap");
        assert_eq!(ov.source_name, "t:0-5000+");
        assert_eq!(ov.sink_name, "q:400-5000+");
        assert_eq!(ov.edge_len, 1400);
    }

    /// Same on the other read: t trimmed at its start, q -> t edge is 1000 + 0 = 1000 + trim.
    #[test]
    fn trimmed_target_start_shifts_query() {
        let reads = [read(0, 5000), read(300, 5000)];
        let ov = classify(&aln(1000, 5000, '+', 0, 4000), &reads).expect("overlap");
        assert_eq!(ov.source_name, "q:0-5000+");
        assert_eq!(ov.edge_len, 1300);
    }

    /// Reverse strand: q[0..4000] aligns to the reverse complement of t[0..4000], so t_rc starts
    /// 1000 bp before q. Trimming q's start by 400 advances the target from its end, giving 1400.
    #[test]
    fn trimmed_query_start_shifts_target_end_on_reverse_strand() {
        let reads = [read(400, 5000), read(0, 5000)];
        let ov = classify(&aln(0, 4000, '-', 0, 4000), &reads).expect("overlap");
        assert_eq!(ov.source_name, "t:0-5000-");
        assert_eq!(ov.edge_len, 1400);
    }

    #[test]
    fn untrimmed_reads_are_unchanged() {
        let reads = [read(0, 5000), read(0, 5000)];
        let ov = classify(&aln(1000, 5000, '+', 0, 4000), &reads).expect("overlap");
        assert_eq!(ov.edge_len, 1000);
    }
}