use std::path::Path;

use trex::illumina::diploid::ParentReferenceParams;
use trex::illumina::mate::MatePairOrientation;
use trex::illumina::multik::MultiKParams;
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
    assert!(out.outputs.fragmentation_path().exists());
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
    let trust = std::fs::read_to_string(out.outputs.trust_path()).unwrap();
    assert!(trust.contains("\"report_only\""));
    assert!(trust.contains("\"trusted_fraction\""));
    assert!(trust.contains("\"correction_candidates\""));
    let simplification = std::fs::read_to_string(out.outputs.simplification_path()).unwrap();
    assert!(simplification.contains("\"scheduler\""));
    assert!(simplification.contains("\"spades_iterative_v1\""));
    assert!(simplification.contains("\"tips\""));
    assert!(simplification.contains("\"diamonds\""));
    let scaffolds = std::fs::read_to_string(out.outputs.scaffolds_path()).unwrap();
    assert!(scaffolds.contains("\"bridge_candidates\""));
    let fragmentation = std::fs::read_to_string(out.outputs.fragmentation_path()).unwrap();
    assert!(fragmentation.contains("\"schema_version\""));
    assert!(fragmentation.contains("\"graph_dead_end_endpoints\""));
    assert_eq!(out.fragmentation_report.summary.contigs, out.contig_count);
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
        multi_k: MultiKParams {
            auto: false,
            ladder: vec![3, 4],
        },
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
    assert!(multi_k.contains("\"dead_end_score\""));
    assert!(multi_k.contains("\"score_terms\""));
}

#[test]
fn multi_k_ladder_uses_selected_k_checkpoint_namespace() {
    let dir = tempfile::tempdir().unwrap();
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let checkpoint_root = dir.path().join("ck");
    let mut params = AssembleParams {
        r1_path: root.join("fixtures/tiny.fq"),
        r2_path: None,
        k: 4,
        trusted_threshold: 1,
        checkpoint_root: Some(checkpoint_root.clone()),
        resume: false,
        strict_checkpoints: false,
        simplify: SimplifyOverrides::default(),
        diploid: DiploidParams::default(),
        multi_k: MultiKParams {
            auto: false,
            ladder: vec![3, 4],
        },
        outputs: outputs_in(dir.path()),
    };

    let first = assemble_illumina(&params).unwrap();
    let selected_dir =
        checkpoint_root.join(format!("selected-k-{}", first.multi_k_selection.selected_k));
    assert!(checkpoint_root.join("preprocess/reads.jsonl").exists());
    assert!(selected_dir.join("counts/kmer_counts.json").exists());
    assert!(selected_dir.join("graph/simplified_dbg.json").exists());
    assert!(selected_dir.join("export/sequences.json").exists());

    params.resume = true;
    let second = assemble_illumina(&params).unwrap();

    assert_eq!(
        first.multi_k_selection.selected_k,
        second.multi_k_selection.selected_k
    );
    assert_eq!(first.unitig_count, second.unitig_count);
    assert_eq!(first.contig_count, second.contig_count);
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
    assert!(second.outputs.stages_path().exists());
    assert!(!second.stage_reports.is_empty());
}

#[test]
fn graph_resume_preserves_simplification_and_mate_sidecars() {
    let dir = tempfile::tempdir().unwrap();
    let r1 = write_tmp(dir.path(), "r1.fq", b"@frag/1\nACGTACGT\n+\nIIIIIIII\n");
    let r2 = write_tmp(dir.path(), "r2.fq", b"@frag/2\nCGTAAAAA\n+\nIIIIIIII\n");
    let ck = dir.path().join("ck");
    let mut params = AssembleParams {
        r1_path: r1,
        r2_path: Some(r2),
        k: 4,
        trusted_threshold: 1,
        checkpoint_root: Some(ck.clone()),
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
    let first = assemble_illumina(&params).unwrap();
    let first_evidence = std::fs::read_to_string(first.outputs.evidence_path()).unwrap();
    let first_simplification =
        std::fs::read_to_string(first.outputs.simplification_path()).unwrap();
    assert!(ck.join("graph/stage_artifacts.json").exists());

    params.resume = true;
    let second = assemble_illumina(&params).unwrap();
    let second_evidence = std::fs::read_to_string(second.outputs.evidence_path()).unwrap();
    let second_simplification =
        std::fs::read_to_string(second.outputs.simplification_path()).unwrap();

    assert_eq!(first.evidence, second.evidence);
    assert_eq!(first.mate_bridge_stats, second.mate_bridge_stats);
    assert_eq!(first.simplify_decisions, second.simplify_decisions);
    assert_eq!(first_evidence, second_evidence);
    assert_eq!(first_simplification, second_simplification);
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
    assert_eq!(stats.report_only_pairs, 0);
    assert_eq!(stats.boosted_edges, 1);
    assert_eq!(stats.candidates.len(), 1);
    assert_eq!(
        stats.candidates[0].orientation,
        MatePairOrientation::R1TailToR2Head
    );
    assert_eq!(
        stats.candidates[0]
            .distance
            .as_ref()
            .map(|distance| distance.estimated_gap_bp),
        Some(-8)
    );
    assert_eq!(
        stats.candidates[0]
            .distance
            .as_ref()
            .map(|distance| distance.confidence),
        Some(50)
    );
    assert_eq!(stats.candidates[0].score, 150);
    assert_eq!(out.scaffold_artifact.bridge_candidates.len(), 1);
    assert_eq!(out.scaffold_artifact.endpoint_join_candidates.len(), 0);
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
    assert!(scaffolds_text.contains("\"orientation\""));
    assert!(scaffolds_text.contains("\"r1_tail_to_r2_head\""));
    assert!(scaffolds_text.contains("\"estimated_gap_bp\""));
    assert!(scaffolds_text.contains("\"constraint_id\""));
    assert!(scaffolds_text.contains("\"kbm000001\""));
    assert!(scaffolds_text.contains("\"from_context\""));
    assert!(scaffolds_text.contains("\"distance_bin\""));
    assert!(scaffolds_text.contains("\"support_histogram\""));
    assert!(scaffolds_text.contains("\"blockers\""));
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
fn mate_endpoint_join_candidates_are_report_only_for_missing_dbg_edges() {
    let with_dir = tempfile::tempdir().unwrap();
    let without_dir = tempfile::tempdir().unwrap();
    let r1 = write_tmp(with_dir.path(), "r1.fq", b"@frag/1\nAAAACGT\n+\nIIIIIII\n");
    let r2 = write_tmp(with_dir.path(), "r2.fq", b"@frag/2\nTTAAGGG\n+\nIIIIIII\n");

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

    let stats = with.mate_bridge_stats.as_ref().expect("mate bridge stats");
    assert_eq!(stats.existing_edge_pairs, 0);
    assert_eq!(stats.report_only_pairs, 1);
    assert_eq!(stats.boosted_edges, 0);
    assert_eq!(
        stats.candidates[0].orientation,
        MatePairOrientation::R1TailToR2Head
    );
    assert_eq!(
        stats.candidates[0]
            .distance
            .as_ref()
            .map(|distance| distance.estimated_gap_bp),
        Some(-6)
    );
    assert_eq!(with.scaffold_artifact.bridge_candidates.len(), 1);
    assert_eq!(with.scaffold_artifact.endpoint_join_candidates.len(), 1);
    assert!(with.scaffold_artifact.paths.is_empty());
    assert!(!with.scaffold_artifact.endpoint_join_candidates[0].accepted);
    assert_eq!(
        with.scaffold_artifact.endpoint_join_candidates[0]
            .rejection_reason
            .as_deref(),
        Some("candidate has fewer than two supporting pairs")
    );
    assert_eq!(
        with.scaffold_artifact.endpoint_join_candidates[0].promotion_stage,
        "report_only_candidate"
    );
    assert_eq!(
        with.scaffold_artifact.endpoint_join_candidates[0].source,
        "mate_pair_endpoint_join_rejected"
    );
    assert_eq!(
        std::fs::read(with.outputs.contigs_path()).unwrap(),
        std::fs::read(without.outputs.contigs_path()).unwrap()
    );
    let scaffolds = std::fs::read_to_string(with.outputs.scaffolds_path()).unwrap();
    assert!(scaffolds.contains("\"promotion_policy\""));
    assert!(scaffolds.contains("\"min_support_pairs\""));
    assert!(scaffolds.contains("\"endpoint_join_candidates\""));
    assert!(scaffolds.contains("\"promotion_stage\""));
    assert!(scaffolds.contains("\"orientation\""));
    assert!(scaffolds.contains("\"estimated_gap_bp\""));
    assert!(scaffolds.contains("\"confidence\""));
    assert!(scaffolds.contains("\"absent_dbg_edge_no_graph_edit\""));
    assert!(scaffolds.contains("\"support_histogram\""));
}

#[test]
fn mate_endpoint_join_acceptance_emits_scaffold_fasta_and_gfa_path() {
    let dir = tempfile::tempdir().unwrap();
    let r1 = write_tmp(
        dir.path(),
        "r1.fq",
        b"@frag1/1\nAAAACGT\n+\nIIIIIII\n@frag2/1\nAAAACGT\n+\nIIIIIII\n",
    );
    let r2 = write_tmp(
        dir.path(),
        "r2.fq",
        b"@frag1/2\nTTAAGGG\n+\nIIIIIII\n@frag2/2\nTTAAGGG\n+\nIIIIIII\n",
    );

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

    let result = assemble_illumina(&params).unwrap();

    assert_eq!(result.scaffold_artifact.endpoint_join_candidates.len(), 1);
    assert!(result.scaffold_artifact.endpoint_join_candidates[0].accepted);
    assert_eq!(
        result.scaffold_artifact.endpoint_join_candidates[0].promotion_stage,
        "scaffold_artifact"
    );
    assert_eq!(
        result.scaffold_artifact.endpoint_join_candidates[0].source,
        "mate_pair_endpoint_join_promoted_sidecar"
    );
    assert_eq!(result.scaffold_artifact.paths.len(), 1);

    let scaffolds_fasta = std::fs::read_to_string(result.outputs.scaffolds_fasta_path()).unwrap();
    assert!(scaffolds_fasta.contains(">scf000001"));
    assert!(scaffolds_fasta.lines().any(|line| !line.starts_with('>')));

    let gfa = std::fs::read_to_string(result.outputs.gfa_path_resolved()).unwrap();
    assert!(gfa.contains("P\tscf000001\t"));
    assert!(gfa.contains("utg000001+\tutg000002-"));
    assert!(gfa.contains("TS:Z:trex-scaffold-sidecar"));
    assert!(gfa.contains("GF:Z:scaffolds.fa"));
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
