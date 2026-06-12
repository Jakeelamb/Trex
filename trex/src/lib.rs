//! Trex: Illumina Phase-1 assembly core (**synchronous** public API).
//!
//! Policy vocabulary lives in repository `CONTEXT.md`. This crate implements
//! ingest, preprocess, *k*-mer counting, de Bruijn graph build and simplification,
//! unitigs/contigs export, and checkpoint hooks.

#![forbid(unsafe_code)]

pub mod dbg;
pub mod error;
pub mod evidence;
pub mod illumina;
pub mod kmer;

pub use error::{CheckpointError, GraphError, IngestError, KmerError, TrexError};
