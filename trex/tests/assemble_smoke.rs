use std::path::Path;

use trex::illumina::diploid::ParentReferenceParams;
use trex::illumina::multik::MultiKParams;
use trex::illumina::pipeline::{
    assemble_illumina, AssembleOutputs, AssembleParams, DiploidParams, SimplifyOverrides,
};
use trex::{IngestError, TrexError};

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
        multi_k: Default::default(),
        outputs: outputs_in(dir.path()),
    };
    let out = assemble_illumina(&params).unwrap();
    assert!(!out.trusted_kmers.is_empty());
    assert_eq!(out.reads.len(), 1);
    assert_eq!(out.simplify_stats.tips_removed, 0);
    assert_eq!(out.simplify_stats.diamond_bubbles_resolved, 0);
    assert!(out.evidence.is_empty());
    assert_eq!(out.mate_bridge_stats, None);
    assert!(out.unitig_count >= 1);
    let uf = out.outputs.unitigs_path();
    assert!(uf.exists());
    assert!(out.outputs.evidence_path().exists());
    assert!(out.outputs.annotations_path().exists());
    assert!(out.outputs.simplification_path().exists());
    assert!(out.outputs.scaffolds_path().exists());
    assert!(out.outputs.audit_json_path().exists());
    assert!(out.outputs.audit_tsv_path().exists());
    assert!(!out.outputs.multi_k_path().exists());
    assert!(!out.multi_k_selection.enabled);
    assert_eq!(out.multi_k_selection.selected_k, params.k);
    assert_eq!(out.audit_report.summary.contigs, out.contig_count);
    assert_eq!(out.graph_annotations.summary.unitig_count, out.unitig_count);
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
        multi_k: Default::default(),
        outputs: outputs_in(dir.path()),
    };
    let out = assemble_illumina(&params).unwrap();
    assert_eq!(out.reads.len(), 2);
    assert!(out.evidence.is_empty());
    assert_eq!(out.mate_bridge_stats, None);
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
        multi_k: Default::default(),
        outputs: outputs_in(dir.path()),
    };
    let out = assemble_illumina(&params).unwrap();
    assert!(out.outputs.unitigs_path().exists());
}

#[test]
fn annotations_sidecar_preserves_phase1_golden_outputs() {
    let dir = tempfile::tempdir().unwrap();
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let params = AssembleParams {
        r1_path: root.join("fixtures/tiny.fq"),
        r2_path: None,
        k: 4,
        trusted_threshold: 1,
        checkpoint_root: None,
        resume: false,
        strict_checkpoints: false,
        simplify: SimplifyOverrides::default(),
        diploid: DiploidParams::default(),
        multi_k: Default::default(),
        outputs: outputs_in(dir.path()),
    };
    let out = assemble_illumina(&params).unwrap();

    for name in ["contigs.fa", "unitigs.fa", "graph.gfa"] {
        let got = std::fs::read(dir.path().join(name)).unwrap();
        let expected = std::fs::read(root.join("fixtures/expected/ref_free_smoke").join(name))
            .unwrap_or_else(|err| panic!("read expected {name}: {err}"));
        assert_eq!(got, expected, "{name} changed after annotation sidecar");
    }
    let annotations = std::fs::read_to_string(out.outputs.annotations_path()).unwrap();
    assert!(annotations.contains("\"summary\""));
    let simplification = std::fs::read_to_string(out.outputs.simplification_path()).unwrap();
    assert!(simplification.contains("\"tips\""));
    assert!(simplification.contains("\"diamonds\""));
    let scaffolds = std::fs::read_to_string(out.outputs.scaffolds_path()).unwrap();
    assert!(scaffolds.contains("\"bridge_candidates\""));
    let audit_json = std::fs::read_to_string(out.outputs.audit_json_path()).unwrap();
    assert!(audit_json.contains("\"low_support_kmers\""));
    let audit_tsv = std::fs::read_to_string(out.outputs.audit_tsv_path()).unwrap();
    assert!(audit_tsv.starts_with("kind\tseverity\tcontig\tstart\tend\tmessage"));
    assert_eq!(out.audit_report.summary.low_support_kmers, 0);
    assert_eq!(out.graph_annotations.summary.node_count, 3);
}

#[test]
fn multi_k_ladder_selects_one_graph_and_writes_sidecar() {
    let dir = tempfile::tempdir().unwrap();
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let params = AssembleParams {
        r1_path: root.join("fixtures/tiny.fq"),
        r2_path: None,
        k: 4,
        trusted_threshold: 1,
        checkpoint_root: None,
        resume: false,
        strict_checkpoints: false,
        simplify: SimplifyOverrides::default(),
        diploid: DiploidParams::default(),
        multi_k: MultiKParams { ladder: vec![3, 4] },
        outputs: outputs_in(dir.path()),
    };
    let out = assemble_illumina(&params).unwrap();

    assert!(out.multi_k_selection.enabled);
    assert_eq!(out.multi_k_selection.requested_k, 4);
    assert!([3, 4].contains(&out.multi_k_selection.selected_k));
    assert_eq!(out.multi_k_selection.candidates.len(), 2);
    assert!(out.outputs.multi_k_path().exists());
    let multi_k = std::fs::read_to_string(out.outputs.multi_k_path()).unwrap();
    assert!(multi_k.contains("\"schema_version\""));
    assert!(multi_k.contains("\"selected_k\""));
}

#[test]
fn multi_k_ladder_rejects_checkpoint_root_until_checkpoint_identity_supports_selected_k() {
    let dir = tempfile::tempdir().unwrap();
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let params = AssembleParams {
        r1_path: root.join("fixtures/tiny.fq"),
        r2_path: None,
        k: 4,
        trusted_threshold: 1,
        checkpoint_root: Some(dir.path().join("ck")),
        resume: false,
        strict_checkpoints: false,
        simplify: SimplifyOverrides::default(),
        diploid: DiploidParams::default(),
        multi_k: MultiKParams { ladder: vec![3, 4] },
        outputs: outputs_in(dir.path()),
    };

    let err = assemble_illumina(&params).expect_err("multi-k checkpoints are not supported yet");
    assert!(matches!(
        err,
        TrexError::Ingest(IngestError::MultiKCheckpointUnsupported)
    ));
}

#[test]
fn annotations_sidecar_regenerates_on_resume() {
    let dir = tempfile::tempdir().unwrap();
    let seq = b"ACGTACGTACGT";
    let qual = vec![b'I'; seq.len()];
    let fq_body = format!(
        "@resume1\n{}\n+\n{}\n",
        std::str::from_utf8(seq).unwrap(),
        std::str::from_utf8(&qual).unwrap()
    );
    let fq = write_tmp(dir.path(), "resume.fq", fq_body.as_bytes());
    let ck = dir.path().join("ck");
    let mut params = AssembleParams {
        r1_path: fq,
        r2_path: None,
        k: 4,
        trusted_threshold: 1,
        checkpoint_root: Some(ck),
        resume: false,
        strict_checkpoints: false,
        simplify: SimplifyOverrides::default(),
        diploid: DiploidParams::default(),
        multi_k: Default::default(),
        outputs: outputs_in(dir.path()),
    };
    let first = assemble_illumina(&params).unwrap();
    let first_annotations = std::fs::read_to_string(first.outputs.annotations_path()).unwrap();

    params.resume = true;
    let second = assemble_illumina(&params).unwrap();
    let second_annotations = std::fs::read_to_string(second.outputs.annotations_path()).unwrap();

    assert_eq!(first_annotations, second_annotations);
    assert_eq!(first.graph_annotations, second.graph_annotations);
}

#[test]
fn phase2_mate_bridge_reports_typed_evidence_counters() {
    let dir = tempfile::tempdir().unwrap();
    let r1 = write_tmp(dir.path(), "r1.fq", b"@frag/1\nACGTACGT\n+\nIIIIIIII\n");
    let r2 = write_tmp(dir.path(), "r2.fq", b"@frag/2\nCGTAAAAA\n+\nIIIIIIII\n");
    let params = AssembleParams {
        r1_path: r1,
        r2_path: Some(r2),
        k: 4,
        trusted_threshold: 1,
        checkpoint_root: None,
        resume: false,
        strict_checkpoints: false,
        simplify: SimplifyOverrides::default(),
        diploid: DiploidParams {
            enabled: true,
            insert_mean_bp: Some(8),
            insert_stddev_bp: None,
            ..Default::default()
        },
        multi_k: Default::default(),
        outputs: outputs_in(dir.path()),
    };
    let out = assemble_illumina(&params).unwrap();
    let stats = out
        .mate_bridge_stats
        .expect("phase2 paired insert prior should run mate bridge");
    assert_eq!(stats.pairs_seen, 1);
    assert_eq!(stats.pairs_with_endpoint_kmers, 1);
    assert_eq!(stats.trusted_endpoint_pairs, 1);
    assert_eq!(stats.existing_edge_pairs, 1);
    assert_eq!(stats.boosted_edges, 1);
    assert_eq!(stats.candidates.len(), 1);
    assert_eq!(out.scaffold_artifact.bridge_candidates.len(), 1);
    assert_eq!(out.evidence.records.len(), 1);
    let record = &out.evidence.records[0];
    assert_eq!(record.support.observed, 1);
    assert_eq!(record.support.eligible, 1);
    assert_eq!(record.support.supporting, 1);
    assert_eq!(record.support.applied, 1);
    assert_eq!(record.counters["boosted_edges"], 1);

    let evidence_text = std::fs::read_to_string(out.outputs.evidence_path()).unwrap();
    assert!(evidence_text.contains("\"mate_bridge_existing_edge\""));
    let scaffolds_text = std::fs::read_to_string(out.outputs.scaffolds_path()).unwrap();
    assert!(scaffolds_text.contains("\"support_pairs\""));
}

#[test]
fn mate_scaffold_sidecar_does_not_change_primary_contigs() {
    let with_dir = tempfile::tempdir().unwrap();
    let without_dir = tempfile::tempdir().unwrap();
    let r1 = write_tmp(
        with_dir.path(),
        "r1.fq",
        b"@frag/1\nACGTACGT\n+\nIIIIIIII\n",
    );
    let r2 = write_tmp(
        with_dir.path(),
        "r2.fq",
        b"@frag/2\nCGTAAAAA\n+\nIIIIIIII\n",
    );

    let with_params = AssembleParams {
        r1_path: r1.clone(),
        r2_path: Some(r2.clone()),
        k: 4,
        trusted_threshold: 1,
        checkpoint_root: None,
        resume: false,
        strict_checkpoints: false,
        simplify: SimplifyOverrides::default(),
        diploid: DiploidParams {
            enabled: true,
            insert_mean_bp: Some(8),
            insert_stddev_bp: None,
            ..Default::default()
        },
        multi_k: Default::default(),
        outputs: outputs_in(with_dir.path()),
    };
    let without_params = AssembleParams {
        r1_path: r1,
        r2_path: Some(r2),
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
            ..Default::default()
        },
        multi_k: Default::default(),
        outputs: outputs_in(without_dir.path()),
    };

    let with = assemble_illumina(&with_params).unwrap();
    let without = assemble_illumina(&without_params).unwrap();

    assert_eq!(with.scaffold_artifact.bridge_candidates.len(), 1);
    assert!(without.scaffold_artifact.bridge_candidates.is_empty());
    assert_eq!(
        std::fs::read(with.outputs.contigs_path()).unwrap(),
        std::fs::read(without.outputs.contigs_path()).unwrap()
    );
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
            ..Default::default()
        },
        multi_k: Default::default(),
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
            ..Default::default()
        },
        multi_k: Default::default(),
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
fn phase2_synthetic_parent_refs_emit_diploid_evidence_without_haplotype_fastas() {
    let dir = tempfile::tempdir().unwrap();
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let fixture = root.join("fixtures/phase2_synthetic/reads.fq");
    let p1 = root.join("fixtures/phase2_synthetic/parent1.fa");
    let p2 = root.join("fixtures/phase2_synthetic/parent2.fa");
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
            parent_references: ParentReferenceParams {
                parent1: Some(p1),
                parent2: Some(p2),
            },
        },
        multi_k: Default::default(),
        outputs: outputs_in(dir.path()),
    };
    let out = assemble_illumina(&params).unwrap();

    assert!(out.outputs.diploid_path().exists());
    assert!(out.diploid_evidence.summary.parent_references_supplied);
    assert!(!out.diploid_evidence.summary.full_haplotype_fasta_claimed);
    assert!(out.diploid_evidence.parent_kmers.is_some());
    assert!(out.diploid_evidence.summary.parent_informative_unitigs > 0);
    let diploid_json = std::fs::read_to_string(out.outputs.diploid_path()).unwrap();
    assert!(diploid_json.contains("\"parent1_only_kmers\""));
    let gfa = std::fs::read_to_string(out.outputs.gfa_path_resolved()).unwrap();
    assert!(gfa.contains("PS:Z:parent-specific-kmer-evidence"));
    assert!(gfa.contains("PS:Z:"));
    assert!(!dir.path().join("hap1.fa").exists());
    assert!(!dir.path().join("hap2.fa").exists());
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
        multi_k: Default::default(),
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
        multi_k: Default::default(),
        outputs: outputs_in(dir.path()),
    };
    let out = assemble_illumina(&params).unwrap();
    assert!(out.unitig_count > 0);
    assert!(
        out.outputs.contigs_path().exists(),
        "assembly should write contig FASTA even when some seed walks are unorientable"
    );
}
