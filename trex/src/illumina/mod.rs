//! Illumina Phase-1 ingest and early pipeline stages.

pub mod checkpoint;
pub mod counts;
pub mod fasta;
pub mod fastq;
pub mod io;
pub mod mate;
pub mod paired;
pub mod phase2_primary;
pub mod pipeline;
pub mod preprocess;
pub mod read;

pub use checkpoint::{CheckpointRoot, GraphCheckpointIdentity};
pub use pipeline::{
    assemble_illumina, AssembleOutputs, AssembleParams, AssembleResult, DiploidParams, SimplifyOverrides,
};
pub use read::Read;
