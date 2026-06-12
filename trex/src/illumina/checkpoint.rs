//! **Phase-1 checkpoint layout anchor**: `preprocess/reads.jsonl`, `counts/kmer_counts.json`,
//! `graph/simplified_dbg.json`, `export/sequences.json` (**unitigs + contigs** after stitching),
//! optional per-stage **`manifest.json`** with **SHA-256** when **strict** mode is on at write time.

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::dbg::graph::DbgGraph;
use crate::error::CheckpointError;
use crate::illumina::read::Read;
use crate::kmer::cmp_dna;

pub type KmerCounts = Vec<(Vec<u8>, u64)>;
pub type SequenceRecords = Vec<(String, Vec<u8>)>;
pub type ExportSequences = (SequenceRecords, SequenceRecords);

#[derive(Debug, Clone)]
pub struct CheckpointRoot {
    pub root: PathBuf,
}

impl CheckpointRoot {
    pub fn new(root: PathBuf) -> Self {
        Self { root }
    }

    pub fn preprocess_dir(&self) -> PathBuf {
        self.root.join("preprocess")
    }

    pub fn preprocess_reads(&self) -> PathBuf {
        self.preprocess_dir().join("reads.jsonl")
    }

    pub fn preprocess_manifest(&self) -> PathBuf {
        self.preprocess_dir().join("manifest.json")
    }

    /// **Phase-2 Illumina** paired layout: **R1** record count before **R2** concatenation (`pair_layout.json`).
    pub fn preprocess_pair_layout(&self) -> PathBuf {
        self.preprocess_dir().join("pair_layout.json")
    }

    pub fn selected_k_root(&self, k: usize) -> Self {
        Self {
            root: self.root.join(format!("selected-k-{k}")),
        }
    }

    pub fn counts_dir(&self) -> PathBuf {
        self.root.join("counts")
    }

    pub fn counts_json(&self) -> PathBuf {
        self.counts_dir().join("kmer_counts.json")
    }

    pub fn counts_manifest(&self) -> PathBuf {
        self.counts_dir().join("manifest.json")
    }

    pub fn graph_dir(&self) -> PathBuf {
        self.root.join("graph")
    }

    pub fn graph_json(&self) -> PathBuf {
        self.graph_dir().join("simplified_dbg.json")
    }

    pub fn graph_manifest(&self) -> PathBuf {
        self.graph_dir().join("manifest.json")
    }

    pub fn export_dir(&self) -> PathBuf {
        self.root.join("export")
    }

    pub fn export_sequences_json(&self) -> PathBuf {
        self.export_dir().join("sequences.json")
    }

    pub fn export_manifest(&self) -> PathBuf {
        self.export_dir().join("manifest.json")
    }
}

/// Remove **export** checkpoint artifacts (called when the simplified **graph** checkpoint is rewritten).
pub fn remove_export_checkpoint(ck: &CheckpointRoot) -> Result<(), CheckpointError> {
    let p = ck.export_dir();
    if p.exists() {
        fs::remove_dir_all(&p)?;
    }
    Ok(())
}

#[derive(Serialize, Deserialize)]
struct Manifest {
    version: u32,
    stage: String,
    sha256: String,
}

#[derive(Serialize, Deserialize)]
struct CountsPayload {
    k: usize,
    entries: Vec<(String, u64)>,
}

#[derive(Serialize, Deserialize)]
struct GraphCheckpointPayload {
    k: usize,
    nodes: Vec<(String, u64)>,
    edges: Vec<(String, String, u64)>,
    #[serde(default)]
    diploid_enabled: bool,
    #[serde(default)]
    diploid_paired_end: bool,
    #[serde(default)]
    diploid_insert_mean_bp: Option<u64>,
    #[serde(default)]
    diploid_insert_stddev_bp: Option<u64>,
    #[serde(default)]
    phase2_mate_bridge_v1: bool,
}

/// Fingerprint for **graph** checkpoint reuse on resume: must match the stored payload or the graph is rebuilt.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct GraphCheckpointIdentity {
    pub diploid_enabled: bool,
    pub diploid_paired_end: bool,
    pub diploid_insert_mean_bp: Option<u64>,
    pub diploid_insert_stddev_bp: Option<u64>,
    /// **Phase-2 Illumina** mate-pair bridge boost applied before simplification (edge weights only).
    pub phase2_mate_bridge_v1: bool,
}

fn decode_acgt_key(label: &str, s: &str) -> Result<Vec<u8>, CheckpointError> {
    for (pos, &byte) in s.as_bytes().iter().enumerate() {
        if !matches!(byte, b'A' | b'C' | b'G' | b'T') {
            return Err(CheckpointError::InvalidGraph(format!(
                "{label}: non-ACGT byte {byte} at position {pos} in `{s}`"
            )));
        }
    }
    Ok(s.as_bytes().to_vec())
}

fn sha256_hex(bytes: &[u8]) -> String {
    let mut h = Sha256::new();
    h.update(bytes);
    format!("{:x}", h.finalize())
}

fn verify_digest(path: &Path, expected_hex: &str) -> Result<(), CheckpointError> {
    let bytes = fs::read(path)?;
    let actual = sha256_hex(&bytes);
    if actual != expected_hex {
        return Err(CheckpointError::DigestMismatch {
            path: path.to_path_buf(),
            expected: expected_hex.to_string(),
            actual,
        });
    }
    Ok(())
}

/// Write preprocess checkpoint; when `strict`, record **SHA-256** of `reads.jsonl` in `manifest.json`.
///
/// When `paired_r1_len` is **`Some`**, writes `pair_layout.json` with **`r1_count`** for resume-time
/// **Phase-2 Illumina** mate bridge identity; omit or pass **`None`** for single-end runs.
pub fn write_preprocess_checkpoint(
    ck: &CheckpointRoot,
    reads: &[Read],
    strict: bool,
    paired_r1_len: Option<usize>,
) -> Result<(), CheckpointError> {
    fs::create_dir_all(ck.preprocess_dir())?;
    let path = ck.preprocess_reads();
    let mut body: Vec<u8> = Vec::new();
    for r in reads {
        let line = serde_json::json!({
            "id": r.id,
            "seq": String::from_utf8_lossy(&r.sequence),
        });
        serde_json::to_writer(&mut body, &line)?;
        body.push(b'\n');
    }
    fs::write(path, &body)?;
    if strict {
        let digest = sha256_hex(&body);
        let manifest = Manifest {
            version: 1,
            stage: "preprocess".into(),
            sha256: digest,
        };
        fs::write(
            ck.preprocess_manifest(),
            serde_json::to_vec_pretty(&manifest)?,
        )?;
    }
    let pair_path = ck.preprocess_pair_layout();
    if let Some(n) = paired_r1_len {
        let payload = serde_json::json!({ "r1_count": n });
        fs::write(&pair_path, serde_json::to_vec_pretty(&payload)?)?;
    } else if pair_path.exists() {
        let _ = fs::remove_file(&pair_path);
    }
    Ok(())
}

/// Load **`pair_layout.json`** if present (**paired-end** preprocess checkpoint).
pub fn load_pair_layout_checkpoint(ck: &CheckpointRoot) -> Result<Option<usize>, CheckpointError> {
    let path = ck.preprocess_pair_layout();
    if !path.exists() {
        return Ok(None);
    }
    #[derive(Deserialize)]
    struct PairLayout {
        r1_count: usize,
    }
    let pl: PairLayout = serde_json::from_slice(&fs::read(&path)?)?;
    Ok(Some(pl.r1_count))
}

/// Load preprocess checkpoint if `reads.jsonl` exists. When `strict`, `manifest.json` must exist and match.
pub fn load_preprocess_checkpoint(
    ck: &CheckpointRoot,
    strict: bool,
) -> Result<Option<Vec<Read>>, CheckpointError> {
    #[derive(Deserialize)]
    struct Record {
        id: String,
        seq: String,
    }

    let path = ck.preprocess_reads();
    if !path.exists() {
        return Ok(None);
    }
    let bytes = fs::read(&path)?;
    let manifest_path = ck.preprocess_manifest();
    if strict {
        if !manifest_path.exists() {
            return Err(CheckpointError::StrictManifestMissing(manifest_path));
        }
        let m: Manifest = serde_json::from_slice(&fs::read(&manifest_path)?)?;
        verify_digest(&path, &m.sha256)?;
    }
    let mut reads = Vec::new();
    for line in bytes.split(|&b| b == b'\n').filter(|l| !l.is_empty()) {
        let rec: Record = serde_json::from_slice(line)?;
        reads.push(Read {
            id: rec.id,
            sequence: rec.seq.into_bytes(),
        });
    }
    Ok(Some(reads))
}

/// Write merged *k*-mer counts (**before** trusted *T* filter) so operators can change *T* on resume.
pub fn write_counts_checkpoint(
    ck: &CheckpointRoot,
    k: usize,
    counts: &[(Vec<u8>, u64)],
    strict: bool,
) -> Result<(), CheckpointError> {
    fs::create_dir_all(ck.counts_dir())?;
    let path = ck.counts_json();
    let entries: Vec<(String, u64)> = counts
        .iter()
        .map(|(key, c)| (String::from_utf8_lossy(key).into_owned(), *c))
        .collect();
    let payload = CountsPayload { k, entries };
    let mut body = serde_json::to_vec(&payload)?;
    body.push(b'\n');
    fs::write(path, &body)?;
    if strict {
        let digest = sha256_hex(&body);
        let manifest = Manifest {
            version: 1,
            stage: "counts".into(),
            sha256: digest,
        };
        fs::write(ck.counts_manifest(), serde_json::to_vec_pretty(&manifest)?)?;
    }
    Ok(())
}

/// Load counts checkpoint if present and **k** matches; otherwise returns **`None`** (caller re-enumerates).
pub fn load_counts_checkpoint(
    ck: &CheckpointRoot,
    strict: bool,
    k: usize,
) -> Result<Option<KmerCounts>, CheckpointError> {
    let path = ck.counts_json();
    if !path.exists() {
        return Ok(None);
    }
    let bytes = fs::read(&path)?;
    let manifest_path = ck.counts_manifest();
    if strict {
        if !manifest_path.exists() {
            return Err(CheckpointError::StrictManifestMissing(manifest_path));
        }
        let m: Manifest = serde_json::from_slice(&fs::read(&manifest_path)?)?;
        verify_digest(&path, &m.sha256)?;
    }
    let payload: CountsPayload = serde_json::from_slice(&bytes)?;
    if payload.k != k {
        return Ok(None);
    }
    let out = payload
        .entries
        .into_iter()
        .map(|(s, c)| (s.into_bytes(), c))
        .collect();
    Ok(Some(out))
}

/// Write simplified **DBG** (**after** tip + bubble simplification, **before** unitigs).
pub fn write_graph_checkpoint(
    ck: &CheckpointRoot,
    graph: &DbgGraph,
    strict: bool,
    identity: &GraphCheckpointIdentity,
) -> Result<(), CheckpointError> {
    remove_export_checkpoint(ck)?;
    fs::create_dir_all(ck.graph_dir())?;
    let mut nodes: Vec<(String, u64)> = graph
        .node_mul
        .iter()
        .map(|(key, mul)| {
            let s = String::from_utf8(key.clone()).map_err(|e| {
                CheckpointError::InvalidGraph(format!("node key is not UTF-8: {e}"))
            })?;
            Ok((s, *mul))
        })
        .collect::<Result<Vec<_>, CheckpointError>>()?;
    nodes.sort_by(|a, b| a.0.cmp(&b.0));

    let mut edges: Vec<(String, String, u64)> = Vec::new();
    for (u, mp) in &graph.adj {
        for (v, w) in mp {
            if cmp_dna(u, v) == Ordering::Less {
                let us = String::from_utf8(u.clone()).map_err(|e| {
                    CheckpointError::InvalidGraph(format!("edge endpoint is not UTF-8: {e}"))
                })?;
                let vs = String::from_utf8(v.clone()).map_err(|e| {
                    CheckpointError::InvalidGraph(format!("edge endpoint is not UTF-8: {e}"))
                })?;
                edges.push((us, vs, *w));
            }
        }
    }
    edges.sort_by(|a, b| (&a.0, &a.1).cmp(&(&b.0, &b.1)));

    let payload = GraphCheckpointPayload {
        k: graph.k,
        nodes,
        edges,
        diploid_enabled: identity.diploid_enabled,
        diploid_paired_end: identity.diploid_paired_end,
        diploid_insert_mean_bp: identity.diploid_insert_mean_bp,
        diploid_insert_stddev_bp: identity.diploid_insert_stddev_bp,
        phase2_mate_bridge_v1: identity.phase2_mate_bridge_v1,
    };
    let mut body = serde_json::to_vec(&payload)?;
    body.push(b'\n');
    let path = ck.graph_json();
    fs::write(path, &body)?;
    if strict {
        let digest = sha256_hex(&body);
        let manifest = Manifest {
            version: 1,
            stage: "graph".into(),
            sha256: digest,
        };
        fs::write(ck.graph_manifest(), serde_json::to_vec_pretty(&manifest)?)?;
    }
    Ok(())
}

/// Load simplified **DBG** checkpoint if present, **`k`** matches, and **diploid resume identity** matches;
/// otherwise **`None`** (caller rebuilds).
pub fn load_graph_checkpoint(
    ck: &CheckpointRoot,
    strict: bool,
    k: usize,
    identity: &GraphCheckpointIdentity,
) -> Result<Option<DbgGraph>, CheckpointError> {
    let path = ck.graph_json();
    if !path.exists() {
        return Ok(None);
    }
    let bytes = fs::read(&path)?;
    let manifest_path = ck.graph_manifest();
    if strict {
        if !manifest_path.exists() {
            return Err(CheckpointError::StrictManifestMissing(manifest_path));
        }
        let m: Manifest = serde_json::from_slice(&fs::read(&manifest_path)?)?;
        verify_digest(&path, &m.sha256)?;
    }
    let payload: GraphCheckpointPayload = serde_json::from_slice(&bytes)?;
    if payload.k != k {
        return Ok(None);
    }
    if payload.diploid_enabled != identity.diploid_enabled
        || payload.diploid_paired_end != identity.diploid_paired_end
        || payload.diploid_insert_mean_bp != identity.diploid_insert_mean_bp
        || payload.diploid_insert_stddev_bp != identity.diploid_insert_stddev_bp
        || payload.phase2_mate_bridge_v1 != identity.phase2_mate_bridge_v1
    {
        return Ok(None);
    }
    let mut node_mul = BTreeMap::new();
    for (s, mul) in payload.nodes {
        let key = decode_acgt_key("node", &s)?;
        node_mul.insert(key, mul);
    }
    let mut edge_rows: Vec<(Vec<u8>, Vec<u8>, u64)> = Vec::with_capacity(payload.edges.len());
    for (us, vs, w) in payload.edges {
        let u = decode_acgt_key("edge.u", &us)?;
        let v = decode_acgt_key("edge.v", &vs)?;
        edge_rows.push((u, v, w));
    }
    let graph = DbgGraph::from_checkpoint_parts(k, node_mul, edge_rows)
        .map_err(|e| CheckpointError::InvalidGraph(e.to_string()))?;
    Ok(Some(graph))
}

#[derive(Serialize, Deserialize)]
struct ExportSeqRecord {
    id: String,
    seq: String,
}

#[derive(Serialize, Deserialize)]
struct ExportSequencesPayload {
    k: usize,
    unitigs: Vec<ExportSeqRecord>,
    contigs: Vec<ExportSeqRecord>,
}

/// Persist stitched **unitigs** and **contigs** (ASCII **ACGT** only) for resume after graph load.
pub fn write_export_checkpoint(
    ck: &CheckpointRoot,
    k: usize,
    unitigs: &[(String, Vec<u8>)],
    contigs: &[(String, Vec<u8>)],
    strict: bool,
) -> Result<(), CheckpointError> {
    fs::create_dir_all(ck.export_dir())?;
    let mut urec = Vec::with_capacity(unitigs.len());
    for (i, (header, seq)) in unitigs.iter().enumerate() {
        let id = if header.is_empty() {
            format!("utg{:06}", i + 1)
        } else {
            header.clone()
        };
        let seqs = String::from_utf8(seq.clone())
            .map_err(|e| CheckpointError::InvalidGraph(format!("export unitig {id}: {e}")))?;
        decode_acgt_key("export.unitig.seq", &seqs)?;
        urec.push(ExportSeqRecord { id, seq: seqs });
    }
    let mut crec = Vec::with_capacity(contigs.len());
    for (i, (header, seq)) in contigs.iter().enumerate() {
        let id = if header.is_empty() {
            format!("ctg{:06}", i + 1)
        } else {
            header.clone()
        };
        let seqs = String::from_utf8(seq.clone())
            .map_err(|e| CheckpointError::InvalidGraph(format!("export contig {id}: {e}")))?;
        decode_acgt_key("export.contig.seq", &seqs)?;
        crec.push(ExportSeqRecord { id, seq: seqs });
    }
    let payload = ExportSequencesPayload {
        k,
        unitigs: urec,
        contigs: crec,
    };
    let mut body = serde_json::to_vec(&payload)?;
    body.push(b'\n');
    let path = ck.export_sequences_json();
    fs::write(path, &body)?;
    if strict {
        let digest = sha256_hex(&body);
        let manifest = Manifest {
            version: 1,
            stage: "export".into(),
            sha256: digest,
        };
        fs::write(ck.export_manifest(), serde_json::to_vec_pretty(&manifest)?)?;
    }
    Ok(())
}

/// Load **export** checkpoint when present and **k** matches.
pub fn load_export_checkpoint(
    ck: &CheckpointRoot,
    strict: bool,
    k: usize,
) -> Result<Option<ExportSequences>, CheckpointError> {
    let path = ck.export_sequences_json();
    if !path.exists() {
        return Ok(None);
    }
    let bytes = fs::read(&path)?;
    let manifest_path = ck.export_manifest();
    if strict {
        if !manifest_path.exists() {
            return Err(CheckpointError::StrictManifestMissing(manifest_path));
        }
        let m: Manifest = serde_json::from_slice(&fs::read(&manifest_path)?)?;
        verify_digest(&path, &m.sha256)?;
    }
    let payload: ExportSequencesPayload = serde_json::from_slice(&bytes)?;
    if payload.k != k {
        return Ok(None);
    }
    let mut unitigs = Vec::with_capacity(payload.unitigs.len());
    for rec in payload.unitigs {
        decode_acgt_key("export.unitig.seq", &rec.seq)?;
        unitigs.push((rec.id, rec.seq.into_bytes()));
    }
    let mut contigs = Vec::with_capacity(payload.contigs.len());
    for rec in payload.contigs {
        decode_acgt_key("export.contig.seq", &rec.seq)?;
        contigs.push((rec.id, rec.seq.into_bytes()));
    }
    Ok(Some((unitigs, contigs)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    use tempfile::tempdir;

    use crate::dbg::graph::DbgGraph;

    #[test]
    fn graph_checkpoint_roundtrip() {
        let dir = tempdir().expect("tempdir");
        let ck = CheckpointRoot::new(dir.path().to_path_buf());
        let k = 3usize;
        let mut nodes = BTreeMap::new();
        for s in ["AAA", "AAC", "ACC"] {
            nodes.insert(s.as_bytes().to_vec(), 2u64);
        }
        let mut g = DbgGraph::new(k, nodes);
        g.add_undirected_edge(b"AAA", b"AAC", 1).expect("edge");
        g.add_undirected_edge(b"AAC", b"ACC", 1).expect("edge");
        write_graph_checkpoint(&ck, &g, false, &GraphCheckpointIdentity::default()).expect("write");
        let g2 = load_graph_checkpoint(&ck, false, k, &GraphCheckpointIdentity::default())
            .expect("load")
            .expect("some");
        assert_eq!(g2.k, g.k);
        assert_eq!(g2.node_mul, g.node_mul);
        assert_eq!(g2.adj, g.adj);
    }

    #[test]
    fn graph_checkpoint_rejects_diploid_identity_mismatch() {
        let dir = tempdir().expect("tempdir");
        let ck = CheckpointRoot::new(dir.path().to_path_buf());
        let k = 3usize;
        let mut nodes = BTreeMap::new();
        for s in ["AAA", "AAC", "ACC"] {
            nodes.insert(s.as_bytes().to_vec(), 2u64);
        }
        let mut g = DbgGraph::new(k, nodes);
        g.add_undirected_edge(b"AAA", b"AAC", 1).expect("edge");
        g.add_undirected_edge(b"AAC", b"ACC", 1).expect("edge");
        let id_saved = GraphCheckpointIdentity {
            diploid_enabled: true,
            diploid_paired_end: true,
            diploid_insert_mean_bp: Some(300),
            diploid_insert_stddev_bp: Some(40),
            phase2_mate_bridge_v1: false,
        };
        write_graph_checkpoint(&ck, &g, false, &id_saved).expect("write");
        let id_other = GraphCheckpointIdentity::default();
        let loaded = load_graph_checkpoint(&ck, false, k, &id_other).expect("load");
        assert!(loaded.is_none());
    }

    #[test]
    fn export_checkpoint_roundtrip() {
        let dir = tempdir().expect("tempdir");
        let ck = CheckpointRoot::new(dir.path().to_path_buf());
        let k = 5usize;
        let unitigs = vec![("utg000001".to_string(), b"ACGTACGTAC".to_vec())];
        let contigs = vec![("ctg000001".to_string(), b"ACGTACGTAC".to_vec())];
        write_export_checkpoint(&ck, k, &unitigs, &contigs, false).expect("write export");
        let (u2, c2) = load_export_checkpoint(&ck, false, k)
            .expect("load")
            .expect("some");
        assert_eq!(u2, unitigs);
        assert_eq!(c2, contigs);
    }
}
