//! Illumina **Phase-1** de Bruijn graph (**DBG**): trusted nodes, overlap edges, simplification, unitigs.

pub mod annotate;
pub mod export;
pub mod graph;
pub mod orient;
pub mod simplify;
pub mod unitig;
pub mod walk;

pub use annotate::{
    annotate_graph, GraphAnnotationSummary, GraphAnnotations, NodeAnnotation, NodeDepthClass,
    UnitigAnnotation,
};
pub use export::{
    contig_path_matches_unitig_primary_path, contig_path_partition_full_unitigs,
    primary_contig_paths_for_gfa, unitig_adjacency_links, write_contigs_fasta, write_gfa1,
    write_unitigs_fasta, GfaWriteOptions, UnitigGfaLink,
};
pub use graph::{build_dbg, DbgGraph};
pub use orient::forward_representatives;
pub use simplify::{
    assert_no_self_loops, plan_tip_clips, remove_diamond_bubbles, remove_diamond_bubbles_ext,
    remove_diamond_bubbles_ext_with_decisions, remove_tips, remove_tips_with_decisions,
    DiploidSimplifyMode, SimplifyDecision, SimplifyDecisionAction, SimplifyDecisionLog,
    SimplifyParams, SimplifyStats,
};
pub use unitig::{extract_unitigs, stitch_sequence};
pub use walk::{reference_contig_paths, ContigWalkTieBreak};
