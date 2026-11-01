//! Repair of stretches that the windows cannot polish.
//!
//! A window is only called if enough read segments cover it. A stretch of the draft that no true
//! read matches (a wrong piece taken from a single layout read, for example) is therefore never
//! corrected, whatever its length: the reads that cross it stop aligning where it starts and
//! continue where it ends. Those reads still show up as two alignments of the same read on the
//! same contig, one on each side. If enough reads bridge the stretch in the same way and clearly
//! outnumber the reads that support the draft, the part of each bridging read between the two
//! sides is cut out and a consensus of these parts replaces the stretch.

use super::*;
use std::collections::HashSet;

/// Number of bridging reads needed to replace a stretch.
const MIN_BRIDGING_READS: usize = 5;
/// The bridging reads must outnumber the reads that support the draft stretch by this factor.
const BRIDGE_RATIO: usize = 3;
/// How far (bp) the alignment ends of bridging reads may differ from the median at each side.
const END_TOLERANCE: usize = 100;
/// Bases of the supported flanks that are replaced together with the stretch, so that the
/// alignment ends (which are fuzzy) lie outside of the replaced part.
const FLANK: usize = 150;
/// Longest stretch (between the two sides) that is repaired.
const MAX_STRETCH: usize = 10_000;

/// An alignment of a read to a contig, in contig coordinates.
#[derive(Clone)]
struct Aln {
    read: String,
    rev: bool,
    qlen: usize,
    /// Query interval in the orientation of the alignment.
    q0: usize,
    /// Start of the alignment in the mapping reference (with flank) and in the contig.
    q1: usize,
    ts: usize,
    te: usize,
    cigar: String,
}

/// A stretch of a contig that is replaced.
pub struct Repair {
    pub contig: usize,
    pub start: usize,
    pub end: usize,
    pub seq: Vec<u8>,
    pub bridging: usize,
    pub supporting_draft: usize,
}

/// Runs of consecutive windows with too few read segments: (contig, first window, last window).
pub fn unsupported_runs(windows: &[Vec<Vec<Seg>>], min_depth: usize) -> Vec<(usize, usize, usize)> {
    let mut runs = Vec::new();
    for (ci, per_contig) in windows.iter().enumerate() {
        let mut start = None;
        for (k, segs) in per_contig.iter().enumerate() {
            match (segs.len() < min_depth, start) {
                (true, None) => start = Some(k),
                (false, Some(s)) => {
                    runs.push((ci, s, k - 1));
                    start = None;
                }
                _ => {}
            }
        }
        if let Some(s) = start {
            runs.push((ci, s, per_contig.len() - 1));
        }
    }
    runs
}

/// Query offset (in the orientation of the alignment) at target position `pos`.
fn query_at(cigar: &str, tstart: usize, q0: usize, pos: usize) -> Option<usize> {
    if pos < tstart {
        return None;
    }
    let (mut t, mut q, mut num) = (tstart, q0, 0usize);
    for ch in cigar.bytes() {
        if ch.is_ascii_digit() {
            num = num * 10 + (ch - b'0') as usize;
            continue;
        }
        let n = num;
        num = 0;
        match ch {
            b'M' | b'=' | b'X' => {
                if pos < t + n {
                    return Some(q + pos - t);
                }
                t += n;
                q += n;
            }
            b'D' | b'N' => {
                if pos < t + n {
                    return Some(q);
                }
                t += n;
            }
            b'I' | b'S' => q += n,
            _ => {}
        }
    }
    (pos == t).then_some(q)
}

fn median(values: &mut [usize]) -> usize {
    values.sort_unstable();
    values[values.len() / 2]
}

/// Consensus of sequences that all cover the same stretch end to end.
fn consensus_of_parts(parts: &[Vec<u8>]) -> Option<Vec<u8>> {
    let mut order: Vec<usize> = (0..parts.len()).collect();
    order.sort_by_key(|&i| parts[i].len());
    let backbone = &parts[order[order.len() / 2]];
    let layers: Vec<Layer> = parts
        .iter()
        .map(|p| Layer { seq: p.clone(), b: 0, e: backbone.len() as u32 })
        .collect();
    poa_window(&mut Poa::new(), &mut DpScratch::default(), backbone, &layers)
}

/// Find the stretches in `runs` that can be repaired and the sequence to replace them with.
/// `paf` is the alignment of the reads to the contigs of this round.
pub fn find_repairs(
    paf: &Path,
    reads_path: &Path,
    layouts: &[Option<Layout>],
    runs: &[(usize, usize, usize)],
) -> Result<Vec<Repair>, Box<dyn Error>> {
    // ---- alignments that end near a run ---------------------------------------------------
    let window_len = |lay: &Layout| lay.bounds[1] - lay.bounds[0];
    let regions: Vec<(usize, usize, usize)> = runs
        .iter()
        .map(|&(ci, k0, k1)| {
            let lay = layouts[ci].as_ref().unwrap();
            (ci, lay.bounds[k0], lay.bounds[k1 + 1])
        })
        .collect();
    let mut near: Vec<Vec<Aln>> = vec![Vec::new(); runs.len()];
    for line in BufReader::new(std::fs::File::open(paf)?).lines() {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            continue;
        }
        let Some(ci) = f[5].strip_prefix('c').and_then(|s| s.parse::<usize>().ok()) else { continue };
        let Some(Some(lay)) = layouts.get(ci) else { continue };
        let num = |i: usize| f[i].parse::<usize>().ok();
        let (Some(qlen), Some(qs), Some(qe), Some(ts), Some(te), Some(matches), Some(block)) =
            (num(1), num(2), num(3), num(7), num(8), num(9), num(10))
        else {
            continue;
        };
        let cigar = f[12..].iter().find_map(|t| t.strip_prefix("cg:Z:"));
        let primary = f[12..].iter().all(|t| t.strip_prefix("tp:A:").is_none_or(|v| v == "P"));
        let Some(cigar) = cigar.filter(|_| primary && block > 0 && matches as f64 / block as f64 >= MIN_IDENTITY) else {
            continue;
        };
        // only alignments to the central copy of the contig; the seam of circular contigs is left alone
        if ts < lay.flank || te > lay.flank + lay.len {
            continue;
        }
        let rev = f[4] == "-";
        let aln = Aln {
            read: f[0].to_string(),
            rev,
            qlen,
            q0: if rev { qlen - qe } else { qs },
            q1: if rev { qlen - qs } else { qe },
            ts: ts - lay.flank,
            te: te - lay.flank,
            cigar: cigar.to_string(),
        };
        for (r, &(rc, start, end)) in regions.iter().enumerate() {
            let w = window_len(layouts[rc].as_ref().unwrap());
            if rc == ci && aln.te + w >= start && aln.ts <= end + w {
                near[r].push(aln.clone());
            }
        }
    }

    // ---- bridging reads per run -----------------------------------------------------------
    struct Candidate {
        contig: usize,
        cut: (usize, usize),
        bridges: Vec<(Aln, Aln)>,
        supporting_draft: usize,
    }
    let mut candidates = Vec::new();
    for (r, alns) in near.iter().enumerate() {
        let (ci, start, end) = regions[r];
        let w = window_len(layouts[ci].as_ref().unwrap());
        let mut by_read: HashMap<&str, Vec<&Aln>> = HashMap::new();
        for a in alns {
            by_read.entry(a.read.as_str()).or_default().push(a);
        }
        let mut pairs: Vec<(&Aln, &Aln)> = Vec::new();
        for list in by_read.values() {
            'read: for a1 in list {
                for a2 in list {
                    // in the orientation of the alignment the read advances along the contig
                    let ordered = a1.rev == a2.rev && a1.q1 <= a2.q0 + END_TOLERANCE;
                    if ordered
                        && a1.te + w >= start
                        && a1.te <= end
                        && a2.ts >= start
                        && a2.ts <= end + w
                        && a1.te <= a2.ts + END_TOLERANCE
                    {
                        pairs.push((a1, a2));
                        break 'read;
                    }
                }
            }
        }
        if pairs.len() < MIN_BRIDGING_READS {
            continue;
        }
        let la = median(&mut pairs.iter().map(|p| p.0.te).collect::<Vec<_>>());
        let ra = median(&mut pairs.iter().map(|p| p.1.ts).collect::<Vec<_>>());
        if ra + END_TOLERANCE < la || ra.saturating_sub(la) > MAX_STRETCH {
            continue;
        }
        let (cl, cr) = (la.saturating_sub(FLANK), ra + FLANK);
        let bridges: Vec<(Aln, Aln)> = pairs
            .into_iter()
            .filter(|(a1, a2)| {
                a1.te.abs_diff(la) <= END_TOLERANCE
                    && a2.ts.abs_diff(ra) <= END_TOLERANCE
                    && a1.ts <= cl
                    && a2.te >= cr
            })
            .map(|(a1, a2)| (a1.clone(), a2.clone()))
            .collect();
        // reads that are aligned through the whole stretch support the draft
        let supporting: HashSet<&str> = alns.iter().filter(|a| a.ts <= la && a.te >= ra).map(|a| a.read.as_str()).collect();
        if bridges.len() < MIN_BRIDGING_READS || bridges.len() < BRIDGE_RATIO * supporting.len().max(1) {
            continue;
        }
        candidates.push(Candidate { contig: ci, cut: (cl, cr), bridges, supporting_draft: supporting.len() });
    }
    if candidates.is_empty() {
        return Ok(Vec::new());
    }

    // ---- sequences of the bridging reads ----------------------------------------------------
    let needed: HashSet<&str> = candidates.iter().flat_map(|c| c.bridges.iter().map(|b| b.0.read.as_str())).collect();
    let mut sequences: HashMap<String, Vec<u8>> = HashMap::new();
    let mut reader = crate::utils::open_fastq_reader(reads_path)?;
    while let Some((header, seq, _, _)) = crate::utils::read_fastq_record(&mut reader)? {
        let name = header[1..].split_whitespace().next().unwrap_or("");
        if needed.contains(name) && !sequences.contains_key(name) {
            sequences.insert(name.to_string(), seq.bytes().map(encode).collect());
        }
    }

    // ---- consensus of the parts between the two sides ------------------------------------
    let mut repairs = Vec::new();
    for c in candidates {
        let (cl, cr) = c.cut;
        let mut parts = Vec::new();
        for (a1, a2) in &c.bridges {
            let Some(read) = sequences.get(&a1.read) else { continue };
            let (Some(qa), Some(qb)) = (query_at(&a1.cigar, a1.ts, a1.q0, cl), query_at(&a2.cigar, a2.ts, a2.q0, cr)) else {
                continue;
            };
            if qb <= qa + 100 || qb > a1.qlen {
                continue;
            }
            let part: Vec<u8> = if a1.rev {
                let (fs, fe) = (a1.qlen - qb, a1.qlen - qa);
                read[fs..fe].iter().rev().map(|&b| if b < 4 { 3 - b } else { 4 }).collect()
            } else {
                read[qa..qb].to_vec()
            };
            parts.push(part);
        }
        if parts.len() < MIN_BRIDGING_READS {
            continue;
        }
        let Some(cons) = consensus_of_parts(&parts) else { continue };
        repairs.push(Repair {
            contig: c.contig,
            start: cl,
            end: cr,
            seq: cons.into_iter().map(decode).collect(),
            bridging: parts.len(),
            supporting_draft: c.supporting_draft,
        });
    }
    Ok(repairs)
}

/// Replace the stretches in the contigs; overlapping repairs of one contig are not applied.
/// Returns the repairs that were applied.
pub fn apply(contigs: &mut [Contig], mut repairs: Vec<Repair>) -> Vec<Repair> {
    repairs.sort_by_key(|r| (r.contig, std::cmp::Reverse(r.start)));
    let mut applied: Vec<Repair> = Vec::new();
    for r in repairs {
        if applied.last().is_some_and(|prev| prev.contig == r.contig && r.end > prev.start) {
            continue;
        }
        contigs[r.contig].seq.splice(r.start..r.end, r.seq.iter().copied());
        applied.push(r);
    }
    applied
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn query_position_follows_the_cigar() {
        // 10 matches, 2 inserted read bases, 5 matches, 3 deleted target bases, 10 matches
        let cigar = "10M2I5M3D10M";
        assert_eq!(query_at(cigar, 100, 7, 100), Some(7));
        assert_eq!(query_at(cigar, 100, 7, 109), Some(16));
        assert_eq!(query_at(cigar, 100, 7, 110), Some(19)); // after the insertion
        assert_eq!(query_at(cigar, 100, 7, 116), Some(24)); // inside the deletion
        assert_eq!(query_at(cigar, 100, 7, 118), Some(24)); // after the deletion
        assert_eq!(query_at(cigar, 100, 7, 119), Some(25));
        assert_eq!(query_at(cigar, 100, 7, 99), None);
        assert_eq!(query_at(cigar, 100, 7, 128), Some(7 + 27));
        assert_eq!(query_at(cigar, 100, 7, 129), None);
    }

    #[test]
    fn runs_of_unsupported_windows() {
        let seg = Seg { read: 0, rev: false, qs: 0, qe: 0, qlen: 0, b: 0, e: 0, ident: 1.0 };
        let w = |n: usize| vec![seg.clone(); n];
        let windows = vec![vec![w(9), w(1), w(0), w(9), w(2)], vec![w(5), w(5)]];
        assert_eq!(unsupported_runs(&windows, 3), vec![(0, 1, 2), (0, 4, 4)]);
    }

    #[test]
    fn consensus_of_noisy_parts_recovers_the_sequence() {
        let mut state = 12345u64;
        let mut next = move || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as usize
        };
        let truth: Vec<u8> = (0..1500).map(|_| (next() % 4) as u8).collect();
        let parts: Vec<Vec<u8>> = (0..9)
            .map(|_| {
                let mut s = Vec::new();
                for &b in &truth {
                    match next() % 100 {
                        0..=2 => s.push((b + 1 + (next() % 3) as u8) % 4),
                        3..=4 => {}
                        5..=6 => {
                            s.push(b);
                            s.push((next() % 4) as u8);
                        }
                        _ => s.push(b),
                    }
                }
                s
            })
            .collect();
        let cons = consensus_of_parts(&parts).expect("consensus");
        let mut prev: Vec<usize> = (0..=truth.len()).collect();
        for (i, a) in cons.iter().enumerate() {
            let mut cur = vec![i + 1];
            for (j, b) in truth.iter().enumerate() {
                cur.push((prev[j] + usize::from(a != b)).min(prev[j + 1] + 1).min(cur[j] + 1));
            }
            prev = cur;
        }
        let diff = prev[truth.len()];
        assert!(diff <= 3, "{diff} differences");
    }

    #[test]
    fn repairs_are_applied_back_to_front_and_overlaps_skipped() {
        let mut contigs = vec![Contig { seq: b"AAAAACCCCCGGGGGTTTTT".to_vec(), circular: false }];
        let rep = |start, end, seq: &[u8]| Repair { contig: 0, start, end, seq: seq.to_vec(), bridging: 5, supporting_draft: 1 };
        let applied = apply(&mut contigs, vec![rep(2, 4, b"XYZ"), rep(10, 12, b"Q"), rep(11, 14, b"W")]);
        assert_eq!(applied.len(), 2);
        // (11..14) is applied first; (10..12) overlaps it and is skipped
        assert_eq!(contigs[0].seq, b"AAXYZACCCCCGWGTTTTT".to_vec());
    }
}