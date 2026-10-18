//! Read-guided contig bridging: a second, contig-level assembly round.
//!
//! All reads are mapped back to the unitigs and unitigs spanned by reads are joined:
//!  * a read only votes for a join if its alignment runs off the end of one contig into
//!    the start of the next (end-anchored evidence);
//!  * a contig end followed by two or more well supported partners is a repeat boundary.
//!    Such contigs are traversed, never used as anchors, and removed once every copy has
//!    been resolved by the reads that span it;
//!  * a bridge is accepted only if it dominates at both of its contig ends;
//!  * the bridge sequence is cut from a real spanning read (the one with the median bridge
//!    length) and replaces the flanks of both contigs, so no gap size has to be guessed.
//!    Contigs that overlap by more than `OVERLAP_JOIN_MIN` bp are joined inside the overlap;
//!  * a linear contig whose two ends are bridged to each other is circularised.
//!
//! The junctions are polished from all reads by the consensus step that follows.

use crate::cli::{ReadType, ScaffoldOpts};
use crate::compress_graph::{CompressedGraph, Unitig, UnitigEdge};
use crate::utils::rev_comp;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Bases of each contig end that are replaced by bridge sequence (minimum).
const FLANK: usize = 500;
/// How far an alignment may stop short of a contig end and still reach it.
const END_SLACK: u32 = 250;
/// Gap estimates further than this (+ a fraction of the gap) from the median are outliers.
const OUTLIER_BASE: i64 = 600;
const OUTLIER_FRAC: f64 = 0.2;
/// Contigs overlapping by more than this are joined inside the overlap (no read sequence is
/// inserted); shorter overlaps and gaps are bridged with a read slice.
const OVERLAP_JOIN_MIN: i64 = 300;
/// Bases trimmed from the first contig's tail when joining inside an overlap.
const OVERLAP_TRIM: usize = 300;
/// Contigs shorter than this are never used as anchors.
const MIN_ANCHOR_LEN: usize = 2 * FLANK + 400;
/// Bases by which a chain must beat the best single alignment of the read.
const MIN_CHAIN_GAIN: i64 = FLANK as i64;
/// Maximum number of alignments per read considered for chaining.
const MAX_HITS_PER_READ: usize = 48;

#[derive(Debug, Default)]
pub struct Stats {
    pub accepted: usize,
    pub circularised: usize,
}

/// An oriented contig (index into the unitig list, reversed). As a "leaving" end it means
/// leaving through the tail of the oriented contig.
type OEnd = (usize, bool);

fn flip((c, r): OEnd) -> OEnd {
    (c, !r)
}

fn sign(rev: bool) -> char {
    if rev { '-' } else { '+' }
}

/// Walk a CIGAR from target position `ts` / query position `q0` and return the coordinate on
/// the other axis reached at coordinate `x` of the chosen axis (`by_target`: x is a target
/// coordinate and a query coordinate is returned, otherwise the reverse).
fn cigar_map(cigar: &str, ts: u32, q0: u32, x: u32, by_target: bool) -> Option<u32> {
    let (mut a, mut b) = if by_target { (ts as u64, q0 as u64) } else { (q0 as u64, ts as u64) };
    let x = x as u64;
    if x == a {
        return Some(b as u32);
    }
    let mut num = 0u64;
    for ch in cigar.bytes() {
        if ch.is_ascii_digit() {
            num = num * 10 + (ch - b'0') as u64;
            continue;
        }
        let n = std::mem::take(&mut num);
        // does the operation consume axis `a` / axis `b`?
        let (ca, cb) = match ch {
            b'M' | b'=' | b'X' => (true, true),
            b'D' | b'N' => (by_target, !by_target),
            b'I' | b'S' => (!by_target, by_target),
            _ => continue,
        };
        if ca && x <= a + n {
            return Some((b + if cb { x - a } else { 0 }) as u32);
        }
        a += if ca { n } else { 0 };
        b += if cb { n } else { 0 };
    }
    None
}

/// One read-to-contig alignment.
#[derive(Clone, Debug)]
struct Hit {
    contig: usize,
    rev: bool,
    qs: u32,
    qe: u32,
    qlen: u32,
    ts: u32,
    te: u32,
    tlen: u32,
    ident: f32,
    cigar: Option<String>,
}

impl Hit {
    fn oriented(&self) -> OEnd {
        (self.contig, self.rev)
    }

    fn len(&self) -> u32 {
        self.qe - self.qs
    }

    /// The read enters the oriented contig through its head.
    fn enter_ok(&self) -> bool {
        let f = FLANK as u32;
        if !self.rev {
            self.ts <= END_SLACK && self.te >= f
        } else {
            self.te + END_SLACK >= self.tlen && self.ts + f <= self.tlen
        }
    }

    /// The read leaves the oriented contig through its tail.
    fn leave_ok(&self) -> bool {
        let f = FLANK as u32;
        if !self.rev {
            self.te + END_SLACK >= self.tlen && self.ts + f <= self.tlen
        } else {
            self.ts <= END_SLACK && self.te >= f
        }
    }

    /// Read position (forward read coordinates) of the oriented contig tail.
    fn leave_pos(&self) -> i64 {
        let tail = if !self.rev { self.tlen - self.te } else { self.ts };
        (self.qe + tail) as i64
    }

    /// Read position (forward read coordinates) of the oriented contig head.
    fn enter_pos(&self) -> i64 {
        let head = if !self.rev { self.ts } else { self.tlen - self.te };
        self.qs as i64 - head as i64
    }

    /// Query start on the strand of the oriented contig.
    fn q0(&self) -> u32 {
        if !self.rev { self.qs } else { self.qlen - self.qe }
    }

    /// Read position (forward read coordinates) aligned to coordinate `x` of the oriented contig.
    fn read_pos_at(&self, x: u32) -> Option<i64> {
        let t = if !self.rev { x } else { self.tlen.checked_sub(x)? };
        if t < self.ts || t > self.te {
            return None;
        }
        let q = cigar_map(self.cigar.as_deref()?, self.ts, self.q0(), t, true)? as i64;
        Some(if !self.rev { q } else { self.qlen as i64 - q })
    }

    /// Oriented-contig coordinate aligned to forward read position `r`.
    fn target_at_read(&self, r: i64) -> Option<u32> {
        let q = if !self.rev { r } else { self.qlen as i64 - r };
        if q < self.q0() as i64 || q > (self.q0() + self.len()) as i64 {
            return None;
        }
        let t = cigar_map(self.cigar.as_deref()?, self.ts, self.q0(), q as u32, false)?;
        if !self.rev { Some(t) } else { self.tlen.checked_sub(t) }
    }

    /// Does this hit continue the read past the tail of `a`? Either it starts at its own
    /// contig head (abutting contigs, gap or repeat in between) or it covers the position
    /// where `a` ends (contigs that overlap in sequence).
    fn continues(&self, a: &Hit) -> bool {
        if self.enter_ok() {
            return true;
        }
        let j = a.leave_pos();
        self.cigar.is_some()
            && (self.qs as i64) < j - FLANK as i64
            && (self.qe as i64) >= j + FLANK as i64
            && self.enter_pos() >= a.enter_pos()
    }
}

fn link_ok(a: &Hit, b: &Hit) -> bool {
    a.leave_ok() && b.continues(a) && !(a.contig == b.contig && a.rev != b.rev) && b.qe > a.qe && b.qs >= a.qs
}

/// Hits of one read that form a run of valid end-to-end links.
struct Chain {
    read: String,
    hits: Vec<Hit>,
}

/// The explanation of a read that covers most of it: the path of linked hits (leaving one
/// contig through its tail, continuing in the next) with the largest union of aligned bases.
/// Every read votes for one path only, so a read crossing a repeat contributes to the copy it
/// actually comes from and not to every copy its repeat sequence aligns to.
fn chain_of_read(read: &str, mut hits: Vec<Hit>) -> Option<Chain> {
    if hits.len() < 2 {
        return None;
    }
    if hits.len() > MAX_HITS_PER_READ {
        hits.sort_by_key(|h| std::cmp::Reverse(h.len()));
        hits.truncate(MAX_HITS_PER_READ);
    }
    hits.sort_by_key(|h| (h.qs, h.qe));
    let n = hits.len();
    // (union length, number of hits, predecessor)
    let mut dp: Vec<(i64, usize, Option<usize>)> = hits.iter().map(|h| (h.len() as i64, 1, None)).collect();
    for i in 0..n {
        for j in 0..i {
            if link_ok(&hits[j], &hits[i]) {
                let add = hits[i].qe as i64 - hits[i].qs.max(hits[j].qe) as i64;
                let cand = (dp[j].0 + add, dp[j].1 + 1, Some(j));
                if (cand.0, cand.1) > (dp[i].0, dp[i].1) {
                    dp[i] = cand;
                }
            }
        }
    }
    let best = (0..n).max_by_key(|&i| (dp[i].0, dp[i].1, std::cmp::Reverse(i)))?;
    // A join is only evidence if no single contig explains the read nearly as well. (Reads
    // inside the overlap of two contigs, or inside a repeat, fit several contigs equally.)
    let best_single = hits.iter().map(|h| h.len() as i64).max()?;
    if dp[best].1 < 2 || dp[best].0 < best_single + MIN_CHAIN_GAIN {
        return None;
    }
    let mut idx = vec![best];
    while let Some(p) = dp[*idx.last()?].2 {
        idx.push(p);
    }
    idx.reverse();
    Some(Chain { read: read.to_string(), hits: idx.into_iter().map(|i| hits[i].clone()).collect() })
}

/// Parse a read-vs-contigs PAF (with `cg:Z:` CIGARs; contigs are named `c<index>`) into chains.
fn parse_chains(paf: &Path, usable: &[bool], opts: &ScaffoldOpts) -> std::io::Result<Vec<Chain>> {
    let mut chains = Vec::new();
    let mut cur = String::new();
    let mut hits: Vec<Hit> = Vec::new();
    for line in BufReader::new(File::open(paf)?).lines() {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            continue;
        }
        if f[0] != cur {
            chains.extend(chain_of_read(&cur, std::mem::take(&mut hits)));
            cur = f[0].to_string();
        }
        let num = |i: usize| f[i].parse::<u32>().ok();
        let (Some(qlen), Some(qs), Some(qe), Some(tlen), Some(ts), Some(te), Some(matches), Some(block)) =
            (num(1), num(2), num(3), num(6), num(7), num(8), num(9), num(10))
        else {
            continue;
        };
        let Some(contig) = f[5].strip_prefix('c').and_then(|s| s.parse::<usize>().ok()) else {
            continue;
        };
        let ident = matches as f64 / block.max(1) as f64;
        if !usable.get(contig).copied().unwrap_or(false)
            || qe.saturating_sub(qs) < opts.min_alignment_len
            || ident < opts.min_identity
        {
            continue;
        }
        hits.push(Hit {
            contig,
            rev: f[4] == "-",
            qs,
            qe,
            qlen,
            ts,
            te,
            tlen,
            ident: ident as f32,
            cigar: f[12..].iter().find_map(|t| t.strip_prefix("cg:Z:")).map(str::to_string),
        });
    }
    chains.extend(chain_of_read(&cur, hits));
    Ok(chains)
}

/// partners[e] = how many reads continue from leaving-end `e` into each oriented contig
/// (counted in both read directions).
fn partner_counts(chains: &[Chain]) -> BTreeMap<OEnd, BTreeMap<OEnd, usize>> {
    let mut m: BTreeMap<OEnd, BTreeMap<OEnd, usize>> = BTreeMap::new();
    for w in chains.iter().flat_map(|c| c.hits.windows(2)) {
        let (a, b) = (w[0].oriented(), w[1].oriented());
        *m.entry(a).or_default().entry(b).or_insert(0) += 1;
        *m.entry(flip(b)).or_default().entry(flip(a)).or_insert(0) += 1;
    }
    m
}

/// Partners with enough support, relative to the best partner.
fn significant_partners(partners: &BTreeMap<OEnd, usize>, opts: &ScaffoldOpts) -> Vec<OEnd> {
    let best = partners.values().copied().max().unwrap_or(0);
    let floor = opts.min_reads.max((opts.dominance * best as f64).ceil() as usize);
    partners.iter().filter(|(_, n)| **n >= floor).map(|(e, _)| *e).collect()
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct BridgeKey {
    a: OEnd,
    b: OEnd,
    /// Repeat contigs traversed between `a` and `b`.
    inner: Vec<OEnd>,
}

impl BridgeKey {
    fn mirror(&self) -> BridgeKey {
        BridgeKey { a: flip(self.b), b: flip(self.a), inner: self.inner.iter().rev().map(|&e| flip(e)).collect() }
    }

    /// Overlapping contigs without anything in between are joined directly.
    fn is_overlap_join(&self, gap: i64) -> bool {
        self.inner.is_empty() && gap < -OVERLAP_JOIN_MIN
    }
}

#[derive(Clone, Debug)]
struct Support {
    read: String,
    /// The read runs in the same direction as the canonical bridge key.
    same_dir: bool,
    a: Hit,
    b: Hit,
    gap: i64,
}

struct Bridge {
    key: BridgeKey,
    supports: Vec<Support>,
    gap: i64,
}

/// Consecutive anchor contigs in a chain (with only repeat contigs in between) form bridge supports.
fn collect_supports(chains: &[Chain], anchors: &[bool]) -> BTreeMap<BridgeKey, Vec<Support>> {
    let mut out: BTreeMap<BridgeKey, Vec<Support>> = BTreeMap::new();
    for ch in chains {
        let idx: Vec<usize> = (0..ch.hits.len()).filter(|&i| anchors[ch.hits[i].contig]).collect();
        let mut seen = BTreeSet::new();
        for w in idx.windows(2) {
            let (i, j) = (w[0], w[1]);
            let fwd = BridgeKey {
                a: ch.hits[i].oriented(),
                b: ch.hits[j].oriented(),
                inner: ch.hits[i + 1..j].iter().map(Hit::oriented).collect(),
            };
            let rev = fwd.mirror();
            if fwd == rev {
                continue; // self-mirrored (fold-back) junction
            }
            let same_dir = fwd < rev;
            let key = if same_dir { fwd } else { rev };
            if !seen.insert(key.clone()) {
                continue; // one vote per read and bridge
            }
            out.entry(key).or_default().push(Support {
                read: ch.read.clone(),
                same_dir,
                a: ch.hits[i].clone(),
                b: ch.hits[j].clone(),
                gap: ch.hits[j].enter_pos() - ch.hits[i].leave_pos(),
            });
        }
    }
    out
}

fn median(v: &mut [i64]) -> i64 {
    v.sort_unstable();
    v[v.len() / 2]
}

/// Number of bases cut from both contig ends around a gap/read-slice bridge.
fn bridge_cut(gap: i64) -> usize {
    FLANK + (-gap).max(0) as usize + 100
}

/// Decide which bridges to apply: enough reads, gap outliers removed, and dominant at both ends.
fn decide_bridges(candidates: BTreeMap<BridgeKey, Vec<Support>>, len: &[usize], opts: &ScaffoldOpts) -> Vec<Bridge> {
    let mut bridges: Vec<Bridge> = Vec::new();
    for (key, mut supports) in candidates {
        let med = median(&mut supports.iter().map(|s| s.gap).collect::<Vec<_>>());
        let tol = OUTLIER_BASE + (OUTLIER_FRAC * med.unsigned_abs() as f64) as i64;
        supports.retain(|s| (s.gap - med).abs() <= tol);
        if supports.len() < opts.min_reads {
            continue;
        }
        let gap = median(&mut supports.iter().map(|s| s.gap).collect::<Vec<_>>());
        let fits = if key.is_overlap_join(gap) {
            len[key.a.0] >= OVERLAP_TRIM + 400 && len[key.b.0] >= (-gap) as usize + 400
        } else {
            let need = 2 * bridge_cut(gap) + 400;
            len[key.a.0] >= need && len[key.b.0] >= need
        };
        if fits {
            bridges.push(Bridge { key, supports, gap });
        }
    }

    // dominance at both contig ends, strongest bridges first
    let mut by_end: BTreeMap<OEnd, Vec<(usize, usize)>> = BTreeMap::new();
    for (i, b) in bridges.iter().enumerate() {
        by_end.entry(b.key.a).or_default().push((b.supports.len(), i));
        by_end.entry(flip(b.key.b)).or_default().push((b.supports.len(), i));
    }
    for list in by_end.values_mut() {
        list.sort_by(|x, y| y.0.cmp(&x.0).then(x.1.cmp(&y.1)));
    }
    let dominant = |end: OEnd, i: usize| {
        let list = &by_end[&end];
        list[0].1 == i && list.get(1).is_none_or(|s| s.0 as f64 <= opts.dominance * list[0].0 as f64)
    };
    let mut order: Vec<usize> = (0..bridges.len()).collect();
    order.sort_by(|&x, &y| bridges[y].supports.len().cmp(&bridges[x].supports.len()).then(bridges[x].key.cmp(&bridges[y].key)));
    let mut used = BTreeSet::new();
    let mut keep = BTreeSet::new();
    for i in order {
        let ends = [bridges[i].key.a, flip(bridges[i].key.b)];
        if ends.iter().all(|&e| dominant(e, i) && !used.contains(&e)) && ends[0] != ends[1] {
            used.extend(ends);
            keep.insert(i);
        }
    }
    bridges.into_iter().enumerate().filter(|(i, _)| keep.contains(i)).map(|(_, b)| b).collect()
}

/// How the sequence of a bridge is obtained.
enum Join {
    /// Join inside the overlap of the two contigs (overlap length in bp).
    Overlap(usize),
    /// Insert `reads[read][x..y]` (reverse complemented if `rc`).
    Slice { read: String, x: usize, y: usize, rc: bool },
}

/// Median overlap length of two contigs from the spanning reads: the read position of `a`'s
/// tail (minus the trim) is mapped into `b` through the CIGAR.
fn overlap_length(b: &Bridge) -> Option<usize> {
    let mut v: Vec<i64> = Vec::new();
    for s in &b.supports {
        let pa = (s.a.tlen as usize).checked_sub(OVERLAP_TRIM)? as u32;
        if let Some(q) = s.a.read_pos_at(pa).and_then(|r| s.b.target_at_read(r)) {
            v.push(OVERLAP_TRIM as i64 + q as i64);
        }
    }
    if v.len() < 2 {
        return None;
    }
    let o = median(&mut v);
    (o > OVERLAP_TRIM as i64).then_some(o as usize)
}

/// The spanning read with the median bridge length, and the slice to cut from it.
fn pick_slice(b: &Bridge) -> Option<Join> {
    let cut = bridge_cut(b.gap);
    let mut valid: Vec<(i64, i64, &Support)> = Vec::new();
    for s in &b.supports {
        let x = s.a.read_pos_at((s.a.tlen as usize).checked_sub(cut)? as u32);
        let y = s.b.read_pos_at(cut as u32);
        if let (Some(x), Some(y)) = (x, y) {
            if y >= x && y - x < 5_000_000 {
                valid.push((x, y, s));
            }
        }
    }
    valid.sort_by(|a, b| (a.1 - a.0).cmp(&(b.1 - b.0)).then(b.2.a.ident.total_cmp(&a.2.a.ident)));
    let &(x, y, s) = valid.get(valid.len() / 2)?;
    Some(Join::Slice { read: s.read.clone(), x: x as usize, y: y as usize, rc: !s.same_dir })
}

/// Sequences of the wanted reads (upper case).
fn load_reads(path: &Path, wanted: &BTreeSet<&str>) -> std::io::Result<HashMap<String, String>> {
    let mut found = HashMap::new();
    let mut reader = crate::utils::open_fastq_reader(path)?;
    while found.len() < wanted.len() {
        let Some((header, seq, _, _)) = crate::utils::read_fastq_record(reader.as_mut())? else { break };
        let name = header[1..].split_whitespace().next().unwrap_or("");
        if wanted.contains(name) {
            found.entry(name.to_string()).or_insert_with(|| seq.to_ascii_uppercase());
        }
    }
    Ok(found)
}

struct AcceptedBridge {
    key: BridgeKey,
    /// Bases removed from the tail of `a` and from the head of `b`.
    tail_cut: usize,
    head_cut: usize,
    /// Sequence inserted between the trimmed contigs (empty when joined inside an overlap).
    slice: String,
}

struct Scaffold {
    seq: String,
    circular: bool,
    /// Source contigs in order.
    parts: Vec<OEnd>,
}

#[derive(Default)]
struct Plan {
    scaffolds: Vec<Scaffold>,
    /// Source contigs that are part of a scaffold, and repeat contigs that are fully resolved.
    removed: BTreeSet<usize>,
    stats: Stats,
}

/// Map the reads to the contigs and decide which contigs to join. The contigs are not modified.
fn plan_bridges(
    seqs: &[String],
    circular: &[bool],
    reads: &Path,
    work_dir: &Path,
    opts: &ScaffoldOpts,
    (threads, read_type): (usize, ReadType),
) -> std::io::Result<Plan> {
    let n = seqs.len();
    let mut plan = Plan::default();
    // only linear contigs of reasonable size can be bridged; circular contigs need no mapping
    let usable: Vec<bool> = (0..n).map(|c| !circular[c] && seqs[c].len() >= 2 * FLANK).collect();
    if !usable.contains(&true) {
        return Ok(plan);
    }

    // ---- map reads to contigs ----
    let ref_path = work_dir.join("bridge_ref.fa");
    let paf_path = work_dir.join("bridge_aln.paf");
    std::fs::write(&ref_path, seqs.iter().enumerate().map(|(i, s)| format!(">c{i}\n{s}\n")).collect::<String>())?;
    let mut args: Vec<&str> = read_type.mapping_args().to_vec();
    args.extend(["-c", "--secondary=no", "-M", "0.9", "-v", "1"]);
    crate::utils::run_minimap2(&args, threads, &ref_path, reads, &paf_path)?;

    // ---- evidence ----
    let chains = parse_chains(&paf_path, &usable, opts)?;
    let _ = std::fs::remove_file(&paf_path);
    let _ = std::fs::remove_file(&ref_path);

    // a contig is an anchor unless one of its ends branches into two or more partners
    let partners = partner_counts(&chains);
    let anchors: Vec<bool> = (0..n)
        .map(|c| {
            usable[c]
                && seqs[c].len() >= MIN_ANCHOR_LEN
                && ![false, true].iter().any(|&r| {
                    partners.get(&(c, r)).is_some_and(|m| significant_partners(m, opts).len() >= 2)
                })
        })
        .collect();
    let len: Vec<usize> = seqs.iter().map(String::len).collect();
    let bridges = decide_bridges(collect_supports(&chains, &anchors), &len, opts);

    // ---- sequences of the accepted bridges ----
    let joins: Vec<Option<Join>> = bridges
        .iter()
        .map(|b| if b.key.is_overlap_join(b.gap) { overlap_length(b).map(Join::Overlap) } else { pick_slice(b) })
        .collect();
    let wanted: BTreeSet<&str> = joins
        .iter()
        .filter_map(|j| if let Some(Join::Slice { read, .. }) = j { Some(read.as_str()) } else { None })
        .collect();
    let read_seqs = load_reads(reads, &wanted)?;
    let accepted: Vec<AcceptedBridge> = bridges
        .into_iter()
        .zip(joins)
        .filter_map(|(b, join)| {
            let (tail_cut, head_cut, slice) = match join? {
                Join::Overlap(o) => (OVERLAP_TRIM, o - OVERLAP_TRIM, String::new()),
                Join::Slice { read, x, y, rc } => {
                    let s = read_seqs.get(&read)?.get(x..y)?;
                    let cut = bridge_cut(b.gap);
                    (cut, cut, if rc { rev_comp(s) } else { s.to_string() })
                }
            };
            Some(AcceptedBridge { key: b.key, tail_cut, head_cut, slice })
        })
        .collect();
    plan.stats.accepted = accepted.len();

    // ---- chain the bridges into paths: next[x] = (entered oriented contig, bridge, mirrored?) ----
    let mut next: BTreeMap<OEnd, (OEnd, usize, bool)> = BTreeMap::new();
    let mut has_pred: BTreeSet<OEnd> = BTreeSet::new();
    for (i, b) in accepted.iter().enumerate() {
        next.insert(b.key.a, (b.key.b, i, false));
        has_pred.insert(b.key.b);
        next.insert(flip(b.key.b), (flip(b.key.a), i, true));
        has_pred.insert(flip(b.key.a));
    }
    // path starts first, then cycles
    let starts: Vec<OEnd> = next.keys().filter(|e| !has_pred.contains(e)).chain(next.keys()).copied().collect();
    let mut visited: BTreeSet<usize> = BTreeSet::new();
    let mut inner_used: BTreeMap<usize, usize> = BTreeMap::new();
    for start in starts {
        if !visited.insert(start.0) {
            continue;
        }
        let (mut nodes, mut links, mut circ, mut cur) = (vec![start], Vec::new(), false, start);
        while let Some(&(nxt, bi, mir)) = next.get(&cur) {
            links.push((bi, mir));
            if nxt == start {
                circ = true;
                break;
            }
            if !visited.insert(nxt.0) {
                links.pop();
                break;
            }
            nodes.push(nxt);
            cur = nxt;
        }
        if links.is_empty() {
            continue;
        }

        // sequence: every contig minus the bases cut at its bridged ends, plus the bridge slices
        let k = nodes.len();
        let mut seq = String::new();
        let mut valid = true;
        for (i, &(c, rev)) in nodes.iter().enumerate() {
            let o = if rev { rev_comp(&seqs[c]) } else { seqs[c].clone() };
            let head = if i > 0 || circ {
                let (bi, mir) = links[(i + links.len() - 1) % links.len()];
                if mir { accepted[bi].tail_cut } else { accepted[bi].head_cut }
            } else {
                0
            };
            let (tail, slice) = if i + 1 < k || circ {
                let (bi, mir) = links[i];
                let b = &accepted[bi];
                if mir { (b.head_cut, rev_comp(&b.slice)) } else { (b.tail_cut, b.slice.clone()) }
            } else {
                (0, String::new())
            };
            if head + tail + 200 > o.len() {
                valid = false;
                break;
            }
            seq.push_str(&o[head..o.len() - tail]);
            seq.push_str(&slice);
        }
        if !valid {
            continue;
        }
        for &(bi, _) in &links {
            for e in &accepted[bi].key.inner {
                *inner_used.entry(e.0).or_insert(0) += 1;
            }
        }
        plan.stats.circularised += circ as usize;
        plan.removed.extend(nodes.iter().map(|e| e.0));
        plan.scaffolds.push(Scaffold { seq, circular: circ, parts: nodes });
    }

    // ---- repeat contigs whose copies are all resolved are dropped ----
    for (&r, &used) in &inner_used {
        if anchors[r] || plan.removed.contains(&r) {
            continue;
        }
        let sig: Vec<Vec<OEnd>> = [false, true]
            .iter()
            .filter_map(|&rev| partners.get(&(r, rev)))
            .map(|m| significant_partners(m, opts))
            .collect();
        let copies = sig.iter().map(Vec::len).max().unwrap_or(0).max(1);
        if sig.iter().flatten().all(|e| anchors[e.0]) && used >= copies {
            plan.removed.insert(r);
        }
    }
    Ok(plan)
}

/// Replace merged/dropped unitigs by the scaffolds and remap the graph edges.
fn apply_plan(graph: &mut CompressedGraph, plan: &Plan) {
    let old = std::mem::take(&mut graph.unitigs);
    let ids: Vec<usize> = old.iter().map(|u| u.id).collect();
    let next_id = ids.iter().max().map_or(0, |m| m + 1);

    // (old unitig id, orientation) -> (new unitig id, orientation) for edge ends
    let mut from_map: BTreeMap<(usize, char), (usize, char)> = BTreeMap::new();
    let mut to_map: BTreeMap<(usize, char), (usize, char)> = BTreeMap::new();
    let mut unitigs: Vec<Unitig> = old.into_iter().enumerate().filter(|(i, _)| !plan.removed.contains(i)).map(|(_, u)| u).collect();
    for (si, s) in plan.scaffolds.iter().enumerate() {
        let id = next_id + si;
        let ((c1, r1), (cn, rn)) = (s.parts[0], s.parts[s.parts.len() - 1]);
        if !s.circular {
            from_map.insert((ids[cn], sign(rn)), (id, '+'));
            from_map.insert((ids[c1], sign(!r1)), (id, '-'));
            to_map.insert((ids[c1], sign(r1)), (id, '+'));
            to_map.insert((ids[cn], sign(!rn)), (id, '-'));
        }
        // The read tiling of the source unitigs no longer describes the scaffold.
        unitigs.push(Unitig { id, members: Vec::new(), fasta_seq: Some(s.seq.clone()), topology: if s.circular { 'c' } else { 'l' } });
    }

    let removed_ids: BTreeSet<usize> = plan.removed.iter().map(|&i| ids[i]).collect();
    let remap = |id: &mut usize, ori: &mut char, map: &BTreeMap<(usize, char), (usize, char)>| -> bool {
        if !removed_ids.contains(id) {
            return true;
        }
        map.get(&(*id, *ori)).map(|&(i, o)| (*id, *ori) = (i, o)).is_some()
    };
    let edges = std::mem::take(&mut graph.edges);
    graph.edges = edges
        .into_iter()
        .filter_map(|mut e: UnitigEdge| {
            (remap(&mut e.from, &mut e.from_ori, &from_map) && remap(&mut e.to, &mut e.to_ori, &to_map)).then_some(e)
        })
        .collect();
    graph.unitigs = unitigs;
}

/// Bridge the unitigs of `graph`. Every round re-maps the reads to the contigs of the previous
/// round, so scaffolds can be extended further.
pub fn scaffold_graph(
    graph: &mut CompressedGraph,
    reads: &Path,
    work_dir: &Path,
    opts: &ScaffoldOpts,
    threads: usize,
    read_type: ReadType,
) -> std::io::Result<Stats> {
    let mut total = Stats::default();
    for round in 1..=opts.rounds {
        if graph.unitigs.is_empty() {
            break;
        }
        let seqs: Vec<String> = graph.unitigs.iter().map(|u| u.fasta_seq.as_deref().unwrap_or("").to_ascii_uppercase()).collect();
        let circular: Vec<bool> = graph.unitigs.iter().map(|u| u.topology == 'c').collect();
        let plan = crate::utils::timed_step(
            &format!("Contig bridging round {round}/{}", opts.rounds),
            || plan_bridges(&seqs, &circular, reads, work_dir, opts, (threads, read_type)),
        )?;
        let scaffold_count = plan.scaffolds.len();
        total.accepted += plan.stats.accepted;
        total.circularised += plan.stats.circularised;
        if plan.scaffolds.is_empty() {
            break;
        }
        apply_plan(graph, &plan);
        println!("Bridging round {round}: {scaffold_count} scaffolds, {} unitigs remain", graph.unitigs.len());
    }
    Ok(total)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn hit(contig: usize, rev: bool, qs: u32, qe: u32, ts: u32, te: u32, tlen: u32) -> Hit {
        Hit { contig, rev, qs, qe, qlen: 100_000, ts, te, tlen, ident: 0.9, cigar: None }
    }

    fn opts() -> ScaffoldOpts {
        ScaffoldOpts { rounds: 1, min_reads: 3, dominance: 0.3, min_alignment_len: 1000, min_identity: 0.8 }
    }

    #[test]
    fn cigar_maps_both_ways() {
        // 10M 2I 10M 3D 10M, target starts at 100, query at 5
        let c = "10M2I10M3D10M";
        assert_eq!(cigar_map(c, 100, 5, 100, true), Some(5));
        assert_eq!(cigar_map(c, 100, 5, 105, true), Some(10));
        assert_eq!(cigar_map(c, 100, 5, 115, true), Some(5 + 12 + 5));
        // inside the deletion the query position stays put
        assert_eq!(cigar_map(c, 100, 5, 122, true), Some(5 + 22));
        assert_eq!(cigar_map(c, 100, 5, 125, true), Some(29));
        assert_eq!(cigar_map(c, 100, 5, 500, true), None);
        // and the other way round: inside the insertion the target position stays put
        assert_eq!(cigar_map(c, 100, 5, 5, false), Some(100));
        assert_eq!(cigar_map(c, 100, 5, 17, false), Some(110));
        assert_eq!(cigar_map(c, 100, 5, 29, false), Some(125));
    }

    #[test]
    fn end_reaching_flags() {
        // forward hit covering the whole tail of a 10 kb contig
        let h = hit(0, false, 0, 4000, 6000, 10_000, 10_000);
        assert!(h.leave_ok() && !h.enter_ok());
        assert_eq!(h.leave_pos(), 4000);
        // reverse hit: the oriented contig is the reverse complement, so a hit at the target
        // start is where the read leaves
        let h = hit(0, true, 0, 4000, 0, 4000, 10_000);
        assert!(h.leave_ok() && !h.enter_ok());
        assert_eq!(h.leave_pos(), 4000);
        // an alignment ending mid-contig reaches neither end
        let h = hit(0, false, 0, 4000, 3000, 7000, 10_000);
        assert!(!h.leave_ok() && !h.enter_ok());
        // touching the end by < FLANK bases is no evidence
        assert!(!hit(0, false, 0, 300, 9_700, 10_000, 10_000).leave_ok());
    }

    #[test]
    fn chains_need_end_to_end_links() {
        let a = hit(0, false, 0, 3000, 7000, 10_000, 10_000);
        let b = hit(1, false, 3100, 7000, 0, 3900, 10_000);
        let chain = chain_of_read("r1", vec![b, a.clone()]).unwrap();
        assert_eq!(chain.hits.len(), 2);
        // second contig only hit mid-contig: no chain
        assert!(chain_of_read("r2", vec![a, hit(1, false, 3100, 7000, 2000, 5900, 10_000)]).is_none());
    }

    fn cand(a: usize, b: usize, n: usize, gap: i64) -> (BridgeKey, Vec<Support>) {
        let h = hit(0, false, 0, 1, 0, 1, 10_000);
        let support = |i: usize, gap: i64| Support { read: format!("r{a}{b}{i}"), same_dir: true, a: h.clone(), b: h.clone(), gap };
        (BridgeKey { a: (a, false), b: (b, false), inner: Vec::new() }, (0..n).map(|i| support(i, gap)).collect())
    }

    #[test]
    fn dominance_and_competition() {
        let lens = vec![20_000; 6];
        let decide = |c: Vec<(BridgeKey, Vec<Support>)>| decide_bridges(c.into_iter().collect(), &lens, &opts());
        // 0->1 (10 reads) competes with 0->2 (3 reads): the weak one loses
        let b = decide(vec![cand(0, 1, 10, 100), cand(0, 2, 3, 100)]);
        assert_eq!(b.len(), 1);
        assert_eq!(b[0].key.b, (1, false));
        // 10 vs 6 reads is ambiguous: neither is applied
        assert!(decide(vec![cand(0, 1, 10, 100), cand(0, 2, 6, 100)]).is_empty());
        // two bridges competing for the *incoming* end of contig 3
        assert!(decide(vec![cand(0, 3, 8, 0), cand(1, 3, 7, 0)]).is_empty());
        // too few reads
        assert!(decide(vec![cand(0, 1, 2, 0)]).is_empty());
    }

    #[test]
    fn outlier_gaps_are_dropped() {
        let (k, mut s) = cand(0, 1, 5, 500);
        s.push(Support { gap: 9_000, ..s[0].clone() });
        let b = decide_bridges([(k, s)].into_iter().collect(), &[20_000; 2], &opts());
        assert_eq!(b.len(), 1);
        assert_eq!(b[0].supports.len(), 5);
        assert_eq!(b[0].gap, 500);
    }

    #[test]
    fn mirror_key_is_involutive() {
        let k = BridgeKey { a: (1, false), b: (4, true), inner: vec![(2, false), (3, true)] };
        assert_eq!(k.mirror().mirror(), k);
        assert_eq!(k.mirror().a, (4, false));
        assert_eq!(k.mirror().inner, vec![(3, false), (2, true)]);
    }

    fn unitig(id: usize) -> Unitig {
        Unitig { id, members: Vec::new(), fasta_seq: Some("ACGT".repeat(10)), topology: 'l' }
    }

    fn uedge(from: usize, from_ori: char, to: usize, to_ori: char) -> UnitigEdge {
        UnitigEdge { from, to, from_ori, to_ori, edge_len: 10, overlap_len: 5, identity: 1.0 }
    }

    #[test]
    fn graph_rewrite_remaps_edges() {
        // 0 -> 1 -> 2, a bridge joins 0+ and 1+: edge 1+ -> 2+ must now leave the scaffold and
        // the former edge 0+ -> 1+ disappears.
        let mut g = CompressedGraph {
            unitigs: vec![unitig(0), unitig(1), unitig(2)],
            edges: vec![uedge(0, '+', 1, '+'), uedge(1, '+', 2, '+'), uedge(2, '-', 1, '-'), uedge(1, '-', 0, '-')],
        };
        let plan = Plan {
            scaffolds: vec![Scaffold { seq: "ACGTACGTACGT".into(), circular: false, parts: vec![(0, false), (1, false)] }],
            removed: [0, 1].into(),
            ..Default::default()
        };
        apply_plan(&mut g, &plan);
        assert_eq!(g.unitigs.iter().map(|u| u.id).collect::<Vec<_>>(), vec![2, 3]);
        let e: Vec<_> = g.edges.iter().map(|e| (e.from, e.from_ori, e.to, e.to_ori)).collect();
        assert_eq!(e, vec![(3, '+', 2, '+'), (2, '-', 3, '-')]);
    }

    /// End-to-end test with minimap2 on a tiny synthetic genome (skipped without minimap2).
    #[test]
    fn bridges_two_contigs_end_to_end() {
        if std::process::Command::new("minimap2").arg("--version").output().is_err() {
            return;
        }
        let mut state = 12345u64;
        let genome: String = (0..30_000)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                ['A', 'C', 'G', 'T'][(state >> 33) as usize & 3]
            })
            .collect();
        let dir = std::env::temp_dir().join(format!("ilesta_bridge_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let reads = dir.join("reads.fq");
        let mut f = File::create(&reads).unwrap();
        for i in 0..20usize {
            let r = &genome[6_000 + i * 150..15_000 + i * 150];
            let r = if i % 2 == 0 { rev_comp(r) } else { r.to_string() };
            write!(f, "@r{i}\n{r}\n+\n{}\n", "5".repeat(r.len())).unwrap();
        }
        // two contigs with a 300 bp gap between them; the second one is reverse complemented
        let seqs = vec![genome[..14_000].to_string(), rev_comp(&genome[14_300..])];
        let plan = plan_bridges(&seqs, &[false; 2], &reads, &dir, &opts(), (1, ReadType::Ont)).unwrap();
        assert_eq!(plan.scaffolds.len(), 1);
        let s = &plan.scaffolds[0];
        assert!(!s.circular);
        assert_eq!(s.parts, vec![(0, false), (1, true)]);
        assert_eq!(s.seq, genome);
        let _ = std::fs::remove_dir_all(&dir);
    }
}