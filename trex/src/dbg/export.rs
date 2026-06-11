//! **GFA 1.0** and **FASTA** export (**Phase-1 FASTA header policy** / **GFA segment naming**).
//! Use path **`-`** for **stdout** (per **Phase-1 export layout**).

use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::dbg::graph::DbgGraph;
use crate::error::GraphError;

fn open_out(path: &Path) -> Result<BufWriter<Box<dyn Write + Send>>, GraphError> {
    if path.as_os_str() == "-" {
        return Ok(BufWriter::new(Box::new(std::io::stdout())));
    }
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let f = File::create(path)?;
    Ok(BufWriter::new(Box::new(f)))
}

/// Write unitigs as FASTA (`utg000001`, …).
pub fn write_unitigs_fasta(path: &Path, sequences: &[(String, Vec<u8>)]) -> Result<(), GraphError> {
    let mut w = open_out(path)?;
    for (i, (header, seq)) in sequences.iter().enumerate() {
        let hid = if header.is_empty() {
            format!("utg{:06}", i + 1)
        } else {
            header.clone()
        };
        writeln!(w, ">{}", hid)?;
        writeln!(w, "{}", String::from_utf8_lossy(seq))?;
    }
    w.flush()?;
    Ok(())
}

/// Write contigs as FASTA (`ctg000001`, …).
pub fn write_contigs_fasta(path: &Path, sequences: &[(String, Vec<u8>)]) -> Result<(), GraphError> {
    let mut w = open_out(path)?;
    for (i, (header, seq)) in sequences.iter().enumerate() {
        let hid = if header.is_empty() {
            format!("ctg{:06}", i + 1)
        } else {
            header.clone()
        };
        writeln!(w, ">{}", hid)?;
        writeln!(w, "{}", String::from_utf8_lossy(seq))?;
    }
    w.flush()?;
    Ok(())
}

/// **Phase-2 Illumina** unitig–unitig adjacency for **`L`** records: forward walk from the **last**
/// *k*-mer of unitig *i* to the **first** *k*-mer of unitig *j* when that edge exists in the simplified graph.
#[derive(Debug, Clone)]
pub struct UnitigGfaLink {
    pub from_utg: usize,
    pub from_orient: char,
    pub to_utg: usize,
    pub to_orient: char,
    pub overlap_cigar: String,
}

/// Build sorted, deduplicated **`L`** link rows (`+` / `+` only in v1) from unitig vertex paths.
pub fn unitig_adjacency_links(graph: &DbgGraph, unitig_paths: &[Vec<Vec<u8>>]) -> Vec<UnitigGfaLink> {
    let k = graph.k;
    if k < 2 || unitig_paths.len() < 2 {
        return Vec::new();
    }
    let overlap_cigar = format!("{}M", k.saturating_sub(1));
    let mut seen = BTreeSet::new();
    let mut out = Vec::new();

    for (i, pi) in unitig_paths.iter().enumerate() {
        if pi.len() < 2 {
            continue;
        }
        let li = pi.last().expect("non-empty");
        for (j, pj) in unitig_paths.iter().enumerate() {
            if i == j || pj.len() < 2 {
                continue;
            }
            let fj = pj.first().expect("non-empty");
            if graph
                .adj
                .get(li)
                .and_then(|m| m.get(fj))
                .copied()
                .unwrap_or(0)
                == 0
            {
                continue;
            }
            let key = (i, '+', j, '+');
            if seen.insert(key) {
                out.push(UnitigGfaLink {
                    from_utg: i + 1,
                    from_orient: '+',
                    to_utg: j + 1,
                    to_orient: '+',
                    overlap_cigar: overlap_cigar.clone(),
                });
            }
        }
    }
    out.sort_by_key(|x| (x.from_utg, x.to_utg));
    out
}

/// When a contig's **vertex path** exactly matches one **unitig** path (same canonical *k*-mers in
/// order, or the reverse order for the opposite strand), return **GFA 1.0 `P`** segment steps
/// (`utg` index is **1-based** to match `S` line names).
pub fn contig_path_matches_unitig_primary_path(
    contig_path: &[Vec<u8>],
    unitig_paths: &[Vec<Vec<u8>>],
) -> Option<Vec<(usize, char)>> {
    for (u, up) in unitig_paths.iter().enumerate() {
        if contig_path == up.as_slice() {
            return Some(vec![(u + 1, '+')]);
        }
        if contig_path.len() == up.len() && !up.is_empty() {
            let reverse_ok = contig_path
                .iter()
                .zip(up.iter().rev())
                .all(|(a, b)| a == b);
            if reverse_ok {
                return Some(vec![(u + 1, '-')]);
            }
        }
    }
    None
}

/// Partition a contig **vertex path** into **full** unitig paths (each step is one complete `utg`
/// forward or reverse), when the contig is a concatenation of whole unitig paths in order.
///
/// Returns **`None`** when no greedy full-unitig cover exists (caller may fall back to exact
/// single-unitig match via [`contig_path_matches_unitig_primary_path`]).
pub fn contig_path_partition_full_unitigs(
    contig_path: &[Vec<u8>],
    unitig_paths: &[Vec<Vec<u8>>],
) -> Option<Vec<(usize, char)>> {
    if contig_path.is_empty() {
        return Some(Vec::new());
    }
    let mut i = 0usize;
    let mut out: Vec<(usize, char)> = Vec::new();
    while i < contig_path.len() {
        let mut best: Option<(usize, usize, char)> = None; // (end_exclusive, utg_0based, orient)
        for (u, up) in unitig_paths.iter().enumerate() {
            if up.is_empty() {
                continue;
            }
            let n = up.len();
            if i + n <= contig_path.len() && contig_path[i..i + n] == up[..] {
                let cand = (i + n, u, '+');
                best = Some(match best {
                    None => cand,
                    Some(b) => pick_longer_unitig_prefix(b, cand),
                });
            }
            if i + n <= contig_path.len() {
                let mut rev_ok = true;
                for t in 0..n {
                    if contig_path[i + t] != up[n - 1 - t] {
                        rev_ok = false;
                        break;
                    }
                }
                if rev_ok {
                    let cand = (i + n, u, '-');
                    best = Some(match best {
                        None => cand,
                        Some(b) => pick_longer_unitig_prefix(b, cand),
                    });
                }
            }
        }
        let (end, u, o) = best?;
        if end == i {
            return None;
        }
        out.push((u + 1, o));
        i = end;
    }
    Some(out)
}

fn pick_longer_unitig_prefix(
    a: (usize, usize, char),
    b: (usize, usize, char),
) -> (usize, usize, char) {
    let la = a.0;
    let lb = b.0;
    match la.cmp(&lb) {
        std::cmp::Ordering::Less => b,
        std::cmp::Ordering::Greater => a,
        std::cmp::Ordering::Equal => {
            if a.1 != b.1 {
                if a.1 < b.1 {
                    a
                } else {
                    b
                }
            } else if a.2 == '+' && b.2 == '-' {
                a
            } else if a.2 == '-' && b.2 == '+' {
                b
            } else {
                a
            }
        }
    }
}

/// Build **`P`**-line payloads (`ctg000001`, …): **multi-unitig** full-path partition when possible,
/// else a single **`utg`** step when the contig path exactly matches one unitig (forward/reverse).
pub fn primary_contig_paths_for_gfa(
    contig_paths: &[Vec<Vec<u8>>],
    unitig_paths: &[Vec<Vec<u8>>],
) -> Vec<(String, Vec<(usize, char)>)> {
    contig_paths
        .iter()
        .enumerate()
        .filter_map(|(i, p)| {
            let segs = contig_path_partition_full_unitigs(p, unitig_paths).or_else(|| {
                contig_path_matches_unitig_primary_path(p, unitig_paths)
            })?;
            if segs.is_empty() {
                None
            } else {
                Some((format!("ctg{:06}", i + 1), segs))
            }
        })
        .collect()
}

/// Minimal **GFA 1.0** with `H` + `S` lines; segment names match FASTA headers.
///
/// When `phase2_illumina_diploid` is true, the header line carries an experimental tag
/// (`XX:Z:trex-phase2-illumina`) so downstream tools can detect **Phase-2 Illumina diploid** exports
/// while parsers that ignore unknown `H`-line tags remain compatible.
///
/// When `diploid_unitig_links` is **`Some`**, append **`L`** lines after **`S`** (unitig graph edges;
/// **Phase-2 Illumina graph summaries** may count these separately from **Phase-1 reference-free metrics**).
///
/// When `primary_contig_paths` is non-empty, append **`P`** lines (**GFA 1.0** primary scaffolded paths
/// over **unitig** `S` segments) after **`L`** lines.
///
/// When **`phase2_unphased_hap_paths`** is true (with **`phase2_illumina_diploid`**), emit a second
/// **`P`** line per primary contig named `p2h000001`, … with the same segment walk and
/// **`TS:Z:trex-unphased-hap-mirror`** (unphased dual-path carrier until richer haplotype walks ship).
pub fn write_gfa1(
    path: &Path,
    segments: &[(String, Vec<u8>)],
    phase2_illumina_diploid: bool,
    diploid_unitig_links: Option<&[UnitigGfaLink]>,
    primary_contig_paths: &[(String, Vec<(usize, char)>)],
    phase2_unphased_hap_paths: bool,
) -> Result<(), GraphError> {
    let mut w = open_out(path)?;
    if phase2_illumina_diploid {
        writeln!(w, "H\tVN:Z:1.0\tXX:Z:trex-phase2-illumina")?;
    } else {
        writeln!(w, "H\tVN:Z:1.0")?;
    }
    for (i, (name, seq)) in segments.iter().enumerate() {
        let sid = if name.is_empty() {
            format!("utg{:06}", i + 1)
        } else {
            name.clone()
        };
        let seqs = String::from_utf8_lossy(seq);
        writeln!(w, "S\t{}\t{}", sid, seqs)?;
    }
    if let Some(links) = diploid_unitig_links {
        for link in links {
            let from_id = format!("utg{:06}", link.from_utg);
            let to_id = format!("utg{:06}", link.to_utg);
            writeln!(
                w,
                "L\t{}\t{}\t{}\t{}\t{}",
                from_id, link.from_orient, to_id, link.to_orient, link.overlap_cigar
            )?;
        }
    }
    for (ctg_name, steps) in primary_contig_paths {
        write!(w, "P\t{}\t", ctg_name)?;
        for (si, (utg_idx, orient)) in steps.iter().enumerate() {
            if si > 0 {
                write!(w, "\t")?;
            }
            write!(w, "utg{:06}{}", utg_idx, orient)?;
        }
        writeln!(w, "\t*")?;
    }
    if phase2_illumina_diploid && phase2_unphased_hap_paths {
        for (ctg_name, steps) in primary_contig_paths {
            let hap_id = ctg_name.replacen("ctg", "p2h", 1);
            write!(w, "P\t{}\t", hap_id)?;
            for (si, (utg_idx, orient)) in steps.iter().enumerate() {
                if si > 0 {
                    write!(w, "\t")?;
                }
                write!(w, "utg{:06}{}", utg_idx, orient)?;
            }
            writeln!(
                w,
                "\t*\tTS:Z:trex-unphased-hap-mirror\tXX:Z:{}",
                ctg_name
            )?;
        }
    }
    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        contig_path_matches_unitig_primary_path, contig_path_partition_full_unitigs,
        primary_contig_paths_for_gfa,
    };

    #[test]
    fn contig_path_unitig_forward_and_reverse() {
        let ut: Vec<Vec<u8>> = vec![b"A".to_vec(), b"B".to_vec()];
        let unitigs: Vec<Vec<Vec<u8>>> = vec![ut.clone()];
        assert_eq!(
            contig_path_matches_unitig_primary_path(&ut, &unitigs),
            Some(vec![(1, '+')])
        );
        let rev: Vec<Vec<u8>> = ut.iter().rev().cloned().collect();
        assert_eq!(
            contig_path_matches_unitig_primary_path(&rev, &unitigs),
            Some(vec![(1, '-')])
        );
    }

    #[test]
    fn partition_two_unitigs_linear() {
        let u0: Vec<Vec<u8>> = vec![b"AA".to_vec(), b"BB".to_vec()];
        let u1: Vec<Vec<u8>> = vec![b"BB".to_vec(), b"CC".to_vec()];
        let u2: Vec<Vec<u8>> = vec![b"CC".to_vec()];
        let unitigs = vec![u0.clone(), u1.clone(), u2.clone()];
        let contig: Vec<Vec<u8>> = vec![b"AA".to_vec(), b"BB".to_vec(), b"CC".to_vec()];
        let p = contig_path_partition_full_unitigs(&contig, &unitigs).expect("partition");
        assert_eq!(p, vec![(1, '+'), (3, '+')]);
    }

    #[test]
    fn primary_paths_names_align_with_contig_order() {
        let ut: Vec<Vec<u8>> = vec![b"x".to_vec(), b"y".to_vec()];
        let unitigs: Vec<Vec<Vec<u8>>> = vec![ut.clone()];
        let c1 = vec![b"a".to_vec()];
        let c2 = ut.clone();
        let out = primary_contig_paths_for_gfa(&[c1, c2], &unitigs);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].0, "ctg000002");
        assert_eq!(out[0].1, vec![(1, '+')]);
    }
}
