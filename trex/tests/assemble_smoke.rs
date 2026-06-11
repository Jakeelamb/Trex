use std::path::Path;

use trex::illumina::pipeline::{
    assemble_illumina, AssembleOutputs, AssembleParams, DiploidParams, SimplifyOverrides,
};

fn write_tmp(dir: &Path, name: &str, content: &[u8]) -> std::path::PathBuf {
    let p = dir.join(name);
    std::fs::write(&p, content).unwrap();
    p
}

fn outputs_in(dir: &Path) -> AssembleOutputs {
    AssembleOutputs {
        out_dir: dir.to_path_buf(),
        unitigs_fasta: Path::new("unitigs.fa").to_path_buf(),
        contigs_fasta: Path::new("contigs.fa").to_path_buf(),
        gfa_path: Path::new("graph.gfa").to_path_buf(),
    }
}

fn load_fasta(path: &Path) -> Vec<Vec<u8>> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut seqs = Vec::new();
    let mut cur = Vec::new();
    for line in text.lines() {
        if line.starts_with('>') {
            if !cur.is_empty() {
                seqs.push(cur);
                cur = Vec::new();
            }
        } else {
            cur.extend(line.trim().as_bytes().iter().map(u8::to_ascii_uppercase));
        }
    }
    if !cur.is_empty() {
        seqs.push(cur);
    }
    seqs
}

fn hamming(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).filter(|(x, y)| x != y).count()
}

fn best_parent_substring_hamming(query: &[u8], parent: &[u8]) -> Option<usize> {
    if query.is_empty() || query.len() > parent.len() {
        return None;
    }
    parent
        .windows(query.len())
        .map(|window| hamming(query, window))
        .min()
}

#[test]
fn assemble_single_end_counts_and_trusted() {
    let dir = tempfile::tempdir().unwrap();
    let seq = b"AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT";
    let qual = vec![b'I'; seq.len()];
    let fq_body = format!(
        "@r1\n{}\n+\n{}\n",
        std::str::from_utf8(seq).unwrap(),
        std::str::from_utf8(&qual).unwrap()
    );
    let fq = write_tmp(dir.path(), "s.fq", fq_body.as_bytes());
    let ck = dir.path().join("ck");
    let params = AssembleParams {
        r1_path: fq,
        r2_path: None,
        k: 4,
        trusted_threshold: 1,
        checkpoint_root: Some(ck),
        resume: false,
        strict_checkpoints: false,
        simplify: SimplifyOverrides::default(),
        diploid: DiploidParams::default(),
        outputs: outputs_in(dir.path()),
    };
    let out = assemble_illumina(&params).unwrap();
    assert!(!out.trusted_kmers.is_empty());
    assert_eq!(out.reads.len(), 1);
    assert!(out.unitig_count >= 1);
    let uf = out.outputs.unitigs_path();
    assert!(uf.exists());
    let fa = std::fs::read_to_string(&uf).unwrap();
    assert!(fa.starts_with(">utg"));
}

#[test]
fn paired_illumina_suffix() {
    let dir = tempfile::tempdir().unwrap();
    let r1 = write_tmp(dir.path(), "r1.fq", b"@frag/1\nACGTACGT\n+\nIIIIIIII\n");
    let r2 = write_tmp(dir.path(), "r2.fq", b"@frag/2\nGCTAGCTA\n+\nIIIIIIII\n");
    let params = AssembleParams {
        r1_path: r1,
        r2_path: Some(r2),
        k: 4,
        trusted_threshold: 1,
        checkpoint_root: None,
        resume: false,
        strict_checkpoints: false,
        simplify: SimplifyOverrides::default(),
        diploid: DiploidParams::default(),
        outputs: outputs_in(dir.path()),
    };
    let out = assemble_illumina(&params).unwrap();
    assert_eq!(out.reads.len(), 2);
}

#[test]
fn assemble_from_fasta() {
    let dir = tempfile::tempdir().unwrap();
    let fa = write_tmp(dir.path(), "s.fa", b">read1\nACACACAC\n");
    let params = AssembleParams {
        r1_path: fa,
        r2_path: None,
        k: 4,
        trusted_threshold: 1,
        checkpoint_root: None,
        resume: false,
        strict_checkpoints: false,
        simplify: SimplifyOverrides::default(),
        diploid: DiploidParams::default(),
        outputs: outputs_in(dir.path()),
    };
    let out = assemble_illumina(&params).unwrap();
    assert!(out.outputs.unitigs_path().exists());
}

#[test]
fn assemble_diploid_emits_phase2_gfa_header_tag() {
    let dir = tempfile::tempdir().unwrap();
    let fq_body = b"@r1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n";
    let fq = write_tmp(dir.path(), "s.fq", fq_body);
    let params = AssembleParams {
        r1_path: fq,
        r2_path: None,
        k: 4,
        trusted_threshold: 1,
        checkpoint_root: None,
        resume: false,
        strict_checkpoints: false,
        simplify: SimplifyOverrides::default(),
        diploid: DiploidParams {
            enabled: true,
            insert_mean_bp: None,
            insert_stddev_bp: None,
        },
        outputs: outputs_in(dir.path()),
    };
    assemble_illumina(&params).unwrap();
    let gfa_path = dir.path().join("graph.gfa");
    let header = std::fs::read_to_string(&gfa_path).unwrap();
    let first = header.lines().next().expect("gfa line");
    assert!(
        first.contains("XX:Z:trex-phase2-illumina"),
        "expected phase2 tag on H line: {first:?}"
    );
    assert!(
        header.lines().any(|l| l.starts_with("P\t")),
        "expected at least one GFA P line: {header:?}"
    );
}

#[test]
fn phase2_synthetic_primary_contigs_match_a_parent() {
    let dir = tempfile::tempdir().unwrap();
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let fixture = root.join("fixtures/phase2_synthetic/reads.fq");
    let params = AssembleParams {
        r1_path: fixture,
        r2_path: None,
        k: 4,
        trusted_threshold: 1,
        checkpoint_root: None,
        resume: false,
        strict_checkpoints: false,
        simplify: SimplifyOverrides::default(),
        diploid: DiploidParams {
            enabled: true,
            insert_mean_bp: None,
            insert_stddev_bp: None,
        },
        outputs: outputs_in(dir.path()),
    };
    assemble_illumina(&params).unwrap();

    let p1 = load_fasta(&root.join("fixtures/phase2_synthetic/parent1.fa"))
        .into_iter()
        .next()
        .unwrap();
    let p2 = load_fasta(&root.join("fixtures/phase2_synthetic/parent2.fa"))
        .into_iter()
        .next()
        .unwrap();
    let contigs = load_fasta(&dir.path().join("contigs.fa"));
    assert!(!contigs.is_empty(), "expected at least one Phase-2 contig");
    for contig in contigs {
        let d1 = best_parent_substring_hamming(&contig, &p1).unwrap_or(usize::MAX);
        let d2 = best_parent_substring_hamming(&contig, &p2).unwrap_or(usize::MAX);
        assert_eq!(
            d1.min(d2),
            0,
            "Phase-2 primary contig must stay parent-consistent: {}",
            String::from_utf8_lossy(&contig)
        );
    }
}

#[test]
fn sparse_linear_read_survives_default_tip_clipping() {
    let dir = tempfile::tempdir().unwrap();
    let seq = b"ACGTTGCATGTCAGTACGATCGTTA";
    let qual = vec![b'I'; seq.len()];
    let fq_body = format!(
        "@sparse1\n{}\n+\n{}\n",
        std::str::from_utf8(seq).unwrap(),
        std::str::from_utf8(&qual).unwrap()
    );
    let fq = write_tmp(dir.path(), "sparse.fq", fq_body.as_bytes());
    let params = AssembleParams {
        r1_path: fq,
        r2_path: None,
        k: 9,
        trusted_threshold: 1,
        checkpoint_root: None,
        resume: false,
        strict_checkpoints: false,
        simplify: SimplifyOverrides::default(),
        diploid: DiploidParams::default(),
        outputs: outputs_in(dir.path()),
    };
    let out = assemble_illumina(&params).unwrap();
    assert!(out.unitig_count > 0);
    assert!(out.contig_count > 0);
}

#[test]
fn real_repeated_read_does_not_abort_contig_walk() {
    let dir = tempfile::tempdir().unwrap();
    let seq = b"TGGTACTGGAGCAGAAGAGCTTTCGGTAGTAGAGCTTGATGGAGTTGGTACTGGAGCAGAAGAGCTTTCAGTGGTAGAGCTTGATGGAGTTGGTACTGG";
    let qual = vec![b'I'; seq.len()];
    let fq_body = format!(
        "@ERR1308583.210\n{}\n+\n{}\n",
        std::str::from_utf8(seq).unwrap(),
        std::str::from_utf8(&qual).unwrap()
    );
    let fq = write_tmp(dir.path(), "yeast_repeat.fq", fq_body.as_bytes());
    let params = AssembleParams {
        r1_path: fq,
        r2_path: None,
        k: 21,
        trusted_threshold: 1,
        checkpoint_root: None,
        resume: false,
        strict_checkpoints: false,
        simplify: SimplifyOverrides::default(),
        diploid: DiploidParams::default(),
        outputs: outputs_in(dir.path()),
    };
    let out = assemble_illumina(&params).unwrap();
    assert!(out.unitig_count > 0);
    assert!(
        out.outputs.contigs_path().exists(),
        "assembly should write contig FASTA even when some seed walks are unorientable"
    );
}
