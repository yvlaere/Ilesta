//! Circularisation of linear contigs whose two ends overlap.
//!
//! When the overlap graph does not close a cycle (low coverage at the junction, a contained
//! read removed at the wrong place), a circular replicon can be a linear walk that runs past its
//! own start: the last tens of kb repeat the first tens of kb. Such a contig is
//! recognised by aligning its head to its tail; if the head aligns to the tail end-to-end
//! (same strand, both ends reached up to a tolerated overhang), the duplicated copy is cut off
//! and the contig is marked circular. The overhang tolerance covers walks that ran on past the
//! ring into unrelated sequence; the overlap thresholds scale down for small contigs so that
//! plasmids of a few kb can be closed too.

use crate::compress_graph::CompressedGraph;
use std::fs;
use std::io;
use std::path::Path;

/// Smallest overlap that is taken as evidence for a ring on a large contig. Shorter shared ends
/// are typically repeats (IS elements, rRNA operons) at the boundary of a contig, not a closed
/// replicon. Small contigs (plasmids) use a smaller, length-scaled minimum, see `min_overlap`.
const MIN_OVERLAP: usize = 10_000;
/// Absolute floor of the overlap, however small the contig is. Contigs shorter than twice this
/// are never examined.
const MIN_SHORT_OVERLAP: usize = 1_000;
/// For contigs shorter than `MIN_OVERLAP / MIN_OVERLAP_FRACTION`, the overlap must cover at
/// least this fraction of the contig.
const MIN_OVERLAP_FRACTION: f64 = 0.25;
/// Largest piece of each end that is compared.
const MAX_WINDOW: usize = 300_000;
/// How far the alignment may stop short of the outermost bases of the contig, whatever the
/// overlap length.
const END_SLACK: usize = 1_000;
/// The walk that produced a ring often runs on past its own start into unrelated or
/// low-quality sequence (a chimeric end read, a repeat copy), so the overlap alignment can stop
/// short of the contig end. Up to this fraction of the overlap length is tolerated as overhang
/// at either end (on top of `END_SLACK`).
const MAX_OVERHANG_FRACTION: f64 = 0.5;
/// Minimum matches / alignment block of the end-to-end alignment.
const MIN_IDENTITY: f64 = 0.8;
/// Minimum identity when the overlap is below `MIN_OVERLAP`: a short shared end is only
/// convincing if it is a near-exact copy.
const MIN_SHORT_IDENTITY: f64 = 0.9;

/// Smallest overlap accepted as evidence for a ring on a contig of `len` bases.
fn min_overlap(len: usize) -> usize {
    let scaled = (len as f64 * MIN_OVERLAP_FRACTION) as usize;
    MIN_OVERLAP.min(scaled).max(MIN_SHORT_OVERLAP)
}

/// From the PAF of `head` (query) against `tail` (target), where `head` is the first
/// `window` bases and `tail` the last `window` bases of a contig of `len` bases, return the
/// range `[start, end)` of the contig that forms the ring, if its ends overlap.
fn ring_range(paf: &str, len: usize, window: usize) -> Option<(usize, usize)> {
    let mut best: Option<(usize, usize, usize)> = None; // (overlap, ring start, ring end)
    let min_overlap = min_overlap(len);
    for line in paf.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 || f[4] != "+" {
            continue;
        }
        let num = |i: usize| f[i].parse::<usize>().ok();
        let (qs, qe, tlen, ts, te, matches, block) =
            (num(2)?, num(3)?, num(6)?, num(7)?, num(8)?, num(9)?, num(10)?);
        let overlap = qe.saturating_sub(qs);
        let min_identity = if overlap < MIN_OVERLAP { MIN_SHORT_IDENTITY } else { MIN_IDENTITY };
        if block == 0 || (matches as f64) < min_identity * block as f64 {
            continue;
        }
        // the head must start at (about) the start of the contig, the tail must end at (about)
        // its end; a walk that overran the ring may carry an overhang past the overlap
        let slack = END_SLACK.max((overlap as f64 * MAX_OVERHANG_FRACTION) as usize);
        if qs > slack || te + slack < tlen || overlap < min_overlap || overlap * 2 > len {
            continue;
        }
        let ring_end = len - window + ts;
        if ring_end <= qs || ring_end - qs < len / 2 {
            continue;
        }
        if best.is_none_or(|b| overlap > b.0) {
            best = Some((overlap, qs, ring_end));
        }
    }
    best.map(|(_, s, e)| (s, e))
}

/// Circularise the linear contigs of `graph` whose ends overlap. Returns how many were changed.
pub fn circularise_self_overlaps(
    graph: &mut CompressedGraph,
    work_dir: &Path,
    threads: usize,
) -> io::Result<usize> {
    let mut changed = 0;
    for u in graph.unitigs.iter_mut() {
        let Some(seq) = u.fasta_seq.as_ref() else { continue };
        let len = seq.len();
        if u.topology == 'c' || len < 2 * MIN_SHORT_OVERLAP {
            continue;
        }
        let window = (len / 2).min(MAX_WINDOW);
        let head = work_dir.join("circ_head.fa");
        let tail = work_dir.join("circ_tail.fa");
        let paf = work_dir.join("circ.paf");
        fs::write(&head, format!(">head\n{}\n", &seq[..window]))?;
        fs::write(&tail, format!(">tail\n{}\n", &seq[len - window..]))?;
        crate::utils::run_minimap2(
            &["-x", "asm20", "-c", "--secondary=no"],
            threads,
            &tail,
            &head,
            &paf,
        )?;
        let text = fs::read_to_string(&paf)?;
        for p in [&head, &tail, &paf] {
            let _ = fs::remove_file(p);
        }
        if let Some((start, end)) = ring_range(&text, len, window) {
            crate::ui::detail(format!(
                "circularising contig_{} ({}): its ends overlap by {}",
                u.id,
                crate::ui::bp(len as u64),
                crate::ui::bp((len - (end - start)) as u64)
            ));
            u.fasta_seq = Some(seq[start..end].to_string());
            u.topology = 'c';
            // the read tiling no longer describes the trimmed sequence
            u.members = Vec::new();
            changed += 1;
        }
    }
    Ok(changed)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_the_ring_from_an_end_to_end_alignment() {
        // 100 kb contig, window 50 kb: head[0..30000] == tail[20000..50000] (tail starts at 50000)
        let paf = "head\t50000\t0\t30000\t+\ttail\t50000\t20000\t50000\t29000\t30000\t60\n";
        assert_eq!(ring_range(paf, 100_000, 50_000), Some((0, 70_000)));
        // an alignment that stops well short of the end of the tail (overhang of 20 kb on a
        // 30 kb overlap, more than MAX_OVERHANG_FRACTION) is a repeat, not a ring
        let paf = "head\t50000\t0\t30000\t+\ttail\t50000\t0\t30000\t29000\t30000\t60\n";
        assert_eq!(ring_range(paf, 100_000, 50_000), None);
        // opposite strand, short overlap and low identity are rejected
        let rev = "head\t50000\t0\t30000\t-\ttail\t50000\t20000\t50000\t29000\t30000\t60\n";
        let short = "head\t50000\t0\t5000\t+\ttail\t50000\t45000\t50000\t4900\t5000\t60\n";
        let loose = "head\t50000\t0\t30000\t+\ttail\t50000\t20000\t50000\t9000\t30000\t60\n";
        for paf in [rev, short, loose] {
            assert_eq!(ring_range(paf, 100_000, 50_000), None);
        }
    }

    /// A ring whose walk ran past its start and then continued into unrelated sequence: the
    /// tail copy of the overlap is followed by an overhang, so the alignment stops short of the
    /// contig end. Real PAF line from a 371,238 bp plasmid contig (reference plasmid 343,098 bp):
    /// head[0..21580] aligns to tail[157474..180626], 4,993 bp before the end of the tail.
    #[test]
    fn circularises_overrun_with_overhang() {
        let paf = "head\t185619\t0\t21580\t+\ttail\t185619\t157474\t180626\t21310\t23227\t60\n";
        // the ring is everything before the tail copy of the overlap
        assert_eq!(ring_range(paf, 371_238, 185_619), Some((0, 343_093)));
        // the same alignment among the unrelated repeat hits minimap2 reports for that contig
        let paf = "\
head\t185619\t0\t21580\t+\ttail\t185619\t157474\t180626\t21310\t23227\t60
head\t185619\t183733\t184934\t-\ttail\t185619\t163679\t164880\t1201\t1201\t60
head\t185619\t34025\t35084\t+\ttail\t185619\t43920\t44979\t1047\t1059\t60
head\t185619\t94317\t95376\t+\ttail\t185619\t43922\t44981\t1044\t1059\t60
head\t185619\t167691\t169050\t+\ttail\t185619\t164988\t166347\t1272\t1359\t1
head\t185619\t36916\t37684\t+\ttail\t185619\t73120\t73888\t768\t768\t60
head\t185619\t182533\t183034\t+\ttail\t185619\t40194\t40695\t498\t501\t60
";
        assert_eq!(ring_range(paf, 371_238, 185_619), Some((0, 343_093)));
        // the overhang tolerance is bounded: a 15 kb overhang after a 21.6 kb overlap is not
        let paf = "head\t185619\t0\t21580\t+\ttail\t185619\t149000\t170604\t21310\t21580\t60\n";
        assert_eq!(ring_range(paf, 371_238, 185_619), None);
    }

    /// A 4,477 bp plasmid contig whose last 1,982 bp are an exact copy of its first 1,982 bp
    /// (reference plasmid 2,495 bp). Below MIN_OVERLAP, but a near-exact copy that covers a large
    /// part of a small contig.
    #[test]
    fn circularises_small_plasmid() {
        let paf = "head\t2238\t0\t1982\t+\ttail\t2238\t256\t2238\t1982\t1982\t60\n";
        assert_eq!(ring_range(paf, 4_477, 2_238), Some((0, 2_495)));
        // a short shared end must be a near-exact copy
        let loose = "head\t2238\t0\t1982\t+\ttail\t2238\t256\t2238\t1600\t1982\t60\n";
        assert_eq!(ring_range(loose, 4_477, 2_238), None);
        // an IS-sized shared end (1.3 kb) on a 20 kb contig is a repeat, not a ring
        let is_end = "head\t10000\t0\t1300\t+\ttail\t10000\t8700\t10000\t1300\t1300\t60\n";
        assert_eq!(ring_range(is_end, 20_000, 10_000), None);
    }

    #[test]
    fn min_overlap_scales_with_contig_length() {
        assert_eq!(min_overlap(1_000_000), MIN_OVERLAP);
        assert_eq!(min_overlap(40_000), MIN_OVERLAP);
        assert_eq!(min_overlap(20_000), 5_000);
        assert_eq!(min_overlap(4_477), 1_119);
        assert_eq!(min_overlap(2_000), MIN_SHORT_OVERLAP);
    }

    /// A 60 kb ring that was walked 20 kb too far (the last 20 kb repeat the first 20 kb with
    /// ~3% noise) is cut back to 60 kb and marked circular; a contig without overlapping ends
    /// is left alone.
    #[test]
    #[ignore = "requires minimap2 on PATH"]
    fn trims_and_marks_an_overrun_contig() {
        let mut state = 99u64;
        let mut next = move || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as usize
        };
        let ring: Vec<u8> = (0..60_000).map(|_| b"ACGT"[next() % 4]).collect();
        let mut dup = Vec::new();
        for &b in &ring[..20_000] {
            match next() % 100 {
                0 => {}
                1 => {
                    dup.push(b);
                    dup.push(b"ACGT"[next() % 4]);
                }
                2 => dup.push(b"ACGT"[next() % 4]),
                _ => dup.push(b),
            }
        }
        let mut overrun = ring.clone();
        overrun.extend_from_slice(&dup);
        let unrelated: Vec<u8> = (0..50_000).map(|_| b"ACGT"[next() % 4]).collect();

        let dir = std::env::temp_dir().join(format!("ilesta-circ-test-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let fasta = dir.join("c.fa");
        let text = format!(
            ">contig_0 topology=l\n{}\n>contig_1 topology=l\n{}\n",
            String::from_utf8_lossy(&overrun),
            String::from_utf8_lossy(&unrelated)
        );
        fs::write(&fasta, text).unwrap();
        let mut g = CompressedGraph::from_fasta(&fasta).unwrap();
        let n = circularise_self_overlaps(&mut g, &dir, 1).unwrap();
        let _ = fs::remove_dir_all(&dir);
        assert_eq!(n, 1);
        assert_eq!(g.unitigs[0].topology, 'c');
        let got = g.unitigs[0].fasta_seq.as_ref().unwrap().len();
        assert!(got.abs_diff(60_000) <= 300, "ring length {got}");
        assert_eq!(g.unitigs[1].topology, 'l');
        assert_eq!(g.unitigs[1].fasta_seq.as_ref().unwrap().len(), 50_000);
    }
}