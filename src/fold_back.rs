// Detection of fold-back reads: nanopore reads whose second half is the reverse complement of the
// first (a sequence followed by its own reverse complement, "inverted duplication" artefacts).
//
// A fold-back read is found by comparing the read with itself: a k-mer of the first arm occurs
// again, reverse-complemented, in the second arm, and all such pairs lie on one anti-diagonal
// (their positions sum to a constant, twice the turning point). Unlike an all-vs-all alignment
// this costs time proportional to the read length only.

/// k-mer length. Two copies of a 13-mer survive the read errors on both arms with a probability
/// of 5-15%, while chance matches inside one 100 kb read are rare (L^2 / 2 / 4^13 ~ 75 per read,
/// spread over every possible turning point).
const K: usize = 13;
/// k-mers that occur more often than this in one read (low complexity, repeats) are ignored.
const MAX_OCCURRENCES: usize = 3;
/// Fewest matching k-mer pairs on the anti-diagonal.
const MIN_PAIRS: usize = 40;
/// At least one matching pair per this many bases of arm.
const MAX_PAIR_SPACING: usize = 100;
/// Shortest arm that counts as a fold-back.
pub const MIN_ARM: usize = 2_000;
/// Largest distance between the two arms (hairpin loop, adapter remnants).
const MAX_GAP: i64 = 500;
/// Half-width of the anti-diagonal band (indels shift the diagonal a little along the arms).
const BAND: i64 = 400;
const BUCKET: i64 = 256;

#[derive(Debug, PartialEq)]
pub struct Hairpin {
    /// Position where the read turns around.
    pub center: usize,
    /// Length of the stretch from the start of the first arm to the end of the second arm.
    pub span: usize,
}

fn base_code(b: u8) -> Option<u32> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Look for a hairpin in `seq`: its two arms are reverse complements of each other and adjacent.
pub fn find_hairpin(seq: &[u8]) -> Option<Hairpin> {
    if seq.len() < 2 * MIN_ARM {
        return None;
    }
    let mask: u32 = (1u32 << (2 * K)) - 1;
    // Forward and reverse-complement k-mers as (code << 32 | position), only the half of all codes
    // selected by a hash (the same codes on both sides, so matching k-mers are still found).
    let mut forward: Vec<u64> = Vec::with_capacity(seq.len() / 2);
    let mut reverse: Vec<u64> = Vec::with_capacity(seq.len() / 2);
    let sampled = |code: u32| code.wrapping_mul(0x9E37_79B1) >> 31 == 0;
    let (mut f, mut r, mut valid) = (0u32, 0u32, 0usize);
    for (p, &b) in seq.iter().enumerate() {
        let Some(c) = base_code(b) else {
            valid = 0;
            continue;
        };
        f = ((f << 2) | c) & mask;
        r = (r >> 2) | ((3 - c) << (2 * (K - 1)));
        valid += 1;
        if valid >= K {
            let start = (p + 1 - K) as u64;
            if sampled(f) {
                forward.push(((f as u64) << 32) | start);
            }
            if sampled(r) {
                reverse.push(((r as u64) << 32) | start);
            }
        }
    }
    forward.sort_unstable();
    reverse.sort_unstable();

    // matching pairs (i, j), i < j, with the anti-diagonal s = i + j + K - 1
    let mut pairs: Vec<(u32, u32)> = Vec::new();
    let (mut a, mut b) = (0, 0);
    while a < forward.len() && b < reverse.len() {
        let (ca, cb) = (forward[a] >> 32, reverse[b] >> 32);
        if ca < cb {
            a += 1;
        } else if cb < ca {
            b += 1;
        } else {
            let a_end = a + forward[a..].iter().take_while(|&&x| x >> 32 == ca).count();
            let b_end = b + reverse[b..].iter().take_while(|&&x| x >> 32 == ca).count();
            if a_end - a <= MAX_OCCURRENCES && b_end - b <= MAX_OCCURRENCES {
                for &x in &forward[a..a_end] {
                    for &y in &reverse[b..b_end] {
                        let (i, j) = (x as u32, y as u32);
                        if (i as usize) + K <= j as usize {
                            pairs.push((i, j));
                        }
                    }
                }
            }
            (a, b) = (a_end, b_end);
        }
    }
    if pairs.len() < MIN_PAIRS {
        return None;
    }

    // the fullest anti-diagonal band
    let s_of = |&(i, j): &(u32, u32)| i as i64 + j as i64 + K as i64 - 1;
    let mut buckets = vec![0usize; (2 * seq.len() as i64 / BUCKET) as usize + 2];
    for p in &pairs {
        buckets[(s_of(p) / BUCKET) as usize] += 1;
    }
    let best = (0..buckets.len() - 1).max_by_key(|&b| buckets[b] + buckets[b + 1])?;
    let mut s_values: Vec<i64> = pairs
        .iter()
        .map(s_of)
        .filter(|&s| s / BUCKET == best as i64 || s / BUCKET == best as i64 + 1)
        .collect();
    s_values.sort_unstable();
    let s_mid = s_values[s_values.len() / 2];
    let band: Vec<&(u32, u32)> = pairs.iter().filter(|p| (s_of(p) - s_mid).abs() <= BAND).collect();
    if band.len() < MIN_PAIRS {
        return None;
    }

    let i_min = band.iter().map(|p| p.0).min()? as i64;
    let i_max = band.iter().map(|p| p.0).max()? as i64 + K as i64;
    let j_min = band.iter().map(|p| p.1).min()? as i64;
    let j_max = band.iter().map(|p| p.1).max()? as i64 + K as i64;
    let arm = (i_max - i_min).min(j_max - j_min);
    if arm < MIN_ARM as i64
        || band.len() * MAX_PAIR_SPACING < arm as usize
        || j_min - i_max > MAX_GAP
    {
        return None;
    }
    Some(Hairpin { center: (s_mid / 2) as usize, span: (j_max - i_min) as usize })
}