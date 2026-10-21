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

    /// A read aligned to its own reverse complement, with the two arms adjacent: the read folds
    /// back on itself (a chimeric hairpin read, second half = reverse complement of the first).
    /// Inverted repeats with a gap between the arms (IS elements, rRNA operons) do not qualify.
    fn is_fold_back(&self) -> bool {
        self.is_self_alignment()
            && self.strand == '-'
            && self.query_end - self.query_start >= MIN_FOLD_BACK_ARM
            && self.target_end - self.target_start >= MIN_FOLD_BACK_ARM
            && self.query_start.max(self.target_start) - self.query_end.min(self.target_end)
                <= MAX_FOLD_BACK_GAP
    }

    /// Position in the read where a fold-back read turns around: the middle of the two arms.
    fn fold_center(&self) -> i64 {
        (self.query_start + self.query_end + self.target_start + self.target_end) / 4
    }

    fn percent_identity(&self) -> f32 {
        (self.num_matching as f32 / self.alignment_block_length as f32) * 100.0
    }
}

/// Shortest self-inverted stretch that marks a read as a fold-back chimera.
const MIN_FOLD_BACK_ARM: i64 = 2_000;
/// Largest distance between the two arms of a fold-back read (hairpin loop, adapter remnants).
const MAX_FOLD_BACK_GAP: i64 = 500;

/// The longest stretch of a read of length `len` between its turning points (0 and `len` are the
/// outer boundaries); the first one if there is a tie.
fn longest_unfolded_stretch(centers: &[i64], len: i64) -> (i64, i64) {
    let mut cuts: Vec<i64> = centers.iter().map(|&c| c.clamp(0, len)).collect();
    cuts.push(0);
    cuts.push(len);
    cuts.sort_unstable();
    cuts.windows(2)
        .map(|w| (w[0], w[1]))
        .max_by_key(|&(a, b)| (b - a, std::cmp::Reverse(a)))
        .unwrap()
}

/// First pass over the PAF: for every fold-back read the stretch to keep, `(start, end)`.
fn fold_back_stretches(paf_in: &std::path::Path) -> std::io::Result<HashMap<String, (i64, i64)>> {
    let mut found: HashMap<String, (i64, Vec<i64>)> = HashMap::new();
    for line in BufReader::new(File::open(paf_in)?).lines() {
        let line = line?;
        let mut f = line.splitn(7, '\t');
        let (Some(q), _, _, _, _, Some(t)) = (f.next(), f.next(), f.next(), f.next(), f.next(), f.next()) else {
            continue;
        };
        if q != t {
            continue;
        }
        if let Some(a) = Alignment::from_line(&line).filter(Alignment::is_fold_back) {
            let entry = found.entry(a.query_name.clone()).or_insert((a.query_length as i64, Vec::new()));
            entry.1.push(a.fold_center());
        }
    }
    Ok(found
        .into_iter()
        .map(|(name, (len, centers))| (name, longest_unfolded_stretch(&centers, len)))
        .collect())
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
    reads: &Vec<Read>,
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
    paf_in: &std::path::Path,
    opts: &crate::cli::OverlapOpts,
) -> Result<AlignmentFilteringOutput, Box<dyn std::error::Error>> {
    let crate::cli::OverlapOpts { min_overlap_length, min_overlap_count, min_percent_identity, overhang_ratio } = opts;
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
    // reads that fold back on themselves (self-alignment on the opposite strand), and the part of
    // each that is kept. Known before the alignments are read: of two alignments between the same
    // reads only the longest is stored, and it must not be the one on the discarded half.
    let kept_stretch = fold_back_stretches(paf_in)?;
    let mut fold_back_skipped: usize = 0;
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

            // alignments that fall in the discarded part of a fold-back read are repeats of the kept part
            let outside_kept = |name: &str, start: i64, end: i64| {
                kept_stretch
                    .get(name)
                    .is_some_and(|&(ks, ke)| end.min(ke) - start.max(ks) < i64::from(*min_overlap_length))
            };
            if outside_kept(&record.query_name, record.query_start, record.query_end)
                || outside_kept(&record.target_name, record.target_start, record.target_end)
            {
                fold_back_skipped += 1;
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
                .for_each(|c| *c = c.saturating_add(1));
            reads[target_id].per_base_coverage[tstart..tend]
                .iter_mut()
                .for_each(|c| *c = c.saturating_add(1));
        }
    }

    // all alignments have been read
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

    // A fold-back read is a chimera of a sequence and its reverse complement: everything after
    // the turning point repeats what came before. Such reads are among the longest ones and end
    // up as the terminal read of a contig, where they replace the real contig end by a mirror
    // image. Keep the longest stretch between turning points (so a long genuine part of the
    // read is not lost) and let the usual clipping of alignments to the trimmed read do the rest.
    let mut fold_back_trimmed = 0usize;
    let mut fold_back_bases = 0u64;
    for (name, &(seg_start, seg_end)) in &kept_stretch {
        let Some(&id) = read_name2read_id.get(name) else { continue };
        let read = &mut reads[id];
        fold_back_bases += read.per_base_coverage.len() as u64 - (seg_end - seg_start) as u64;
        let start = read.coverage_start.max(seg_start as u32);
        let end = read.coverage_end.min(seg_end as u32);
        if end > start {
            read.coverage_start = start;
            read.coverage_end = end;
        } else {
            read.coverage_start = 0;
            read.coverage_end = 0;
        }
        fold_back_trimmed += 1;
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
    println!(
        "Alignment filtering: {} reads, {} alignments retained; {} self-alignments, {} short alignments, {} low-identity alignments skipped",
        reads.len(),
        alignments.len(),
        self_alignments_skipped,
        alignment_length_skipped,
        percent_identity_skipped
    );
    println!(
        "Fold-back (self-inverted) reads trimmed: {} ({} bases removed, {} alignments to the discarded part skipped)",
        fold_back_trimmed, fold_back_bases, fold_back_skipped
    );
    println!(
        "Alignments clipped to coverage-trimmed reads: {} ({} bases)",
        CLIPPED.load(std::sync::atomic::Ordering::Relaxed),
        CLIPPED_BASES.load(std::sync::atomic::Ordering::Relaxed)
    );
    println!(
        "Graph input: {} reads, {} overlaps ({} contained, {} low-coverage reads removed)",
        graph_reads.len(),
        overlaps.len(),
        contained_reads.len(),
        low_coverage_reads.len()
    );

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
        classify_alignment(a, 0, 1, &mut overlaps, 0.2, &reads.to_vec(), 500);
        overlaps.remove(&(0, 1))
    }

    impl Clone for Read {
        fn clone(&self) -> Self {
            Read { per_base_coverage: Vec::new(), coverage_start: self.coverage_start, coverage_end: self.coverage_end }
        }
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