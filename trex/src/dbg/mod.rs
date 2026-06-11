//! Illumina **Phase-1** de Bruijn graph (**DBG**): trusted nodes, overlap edges, simplification, unitigs.

pub mod export;
pub mod graph;
pub mod orient;
pub mod simplify;
pub mod unitig;
pub mod walk;

pub use export::{
    contig_path_matches_unitig_primary_path, contig_path_partition_full_unitigs,
    primary_contig_paths_for_gfa, unitig_adjacency_links, write_contigs_fasta, write_gfa1,
    write_unitigs_fasta, UnitigGfaLink,
};
pub use graph::{build_dbg, DbgGraph};
pub use orient::forward_representatives;
pub use simplify::{
    assert_no_self_loops, remove_diamond_bubbles, remove_diamond_bubbles_ext, remove_tips,
    DiploidSimplifyMode, SimplifyParams,
};
pub use unitig::{extract_unitigs, stitch_sequence};
pub use walk::{reference_contig_paths, ContigWalkTieBreak};
