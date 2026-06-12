//! `thiserror`-based errors (**`#[non_exhaustive]`** on public enums).

/// Top-level crate error (aggregate).
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum TrexError {
    #[error(transparent)]
    Ingest(#[from] IngestError),
    #[error(transparent)]
    Kmer(#[from] KmerError),
    #[error(transparent)]
    Checkpoint(#[from] CheckpointError),
    #[error(transparent)]
    Graph(#[from] GraphError),
}

/// de Bruijn graph construction and simplification.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum GraphError {
    #[error("DBG self-loop on canonical k-mer node")]
    SelfLoop,
    #[error("incompatible forward k-mer assignments for canonical node")]
    OrientationConflict,
    #[error("simplified graph still contains a self-adjacency")]
    SimplifiedSelfLoop,
    #[error("I/O error writing assembly output: {0}")]
    Io(#[from] std::io::Error),
    #[error("JSON error writing assembly output: {0}")]
    Json(#[from] serde_json::Error),
}

/// FASTQ ingest, pairing, and preprocess violations.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum IngestError {
    #[error("empty read after preprocess (read id `{0}`)")]
    EmptyRead(String),
    #[error("k ({k}) exceeds shortest post-preprocess read length ({shortest})")]
    KTooLarge { k: usize, shortest: usize },
    #[error("paired read count mismatch: R1 has {r1} reads, R2 has {r2} reads")]
    PairCountMismatch { r1: usize, r2: usize },
    #[error("paired read id mismatch at index {index}: R1=`{r1_id}` R2=`{r2_id}`")]
    PairIdMismatch {
        index: usize,
        r1_id: String,
        r2_id: String,
    },
    #[error(
        "invalid nucleotide byte {byte} (expected ACGTN/IUPAC) in read `{id}` at position {pos}"
    )]
    InvalidNucleotide { id: String, pos: usize, byte: u8 },
    #[error("FASTQ parse error: {0}")]
    FastqFormat(String),
    #[error("`--resume` requires `--checkpoint-root`")]
    ResumeRequiresCheckpointRoot,
    #[error("multi-k selection does not support checkpoints yet; run without `--checkpoint-root`")]
    MultiKCheckpointUnsupported,
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

/// *k*-mer construction and counting.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum KmerError {
    #[error("k must be positive, got {0}")]
    KZero(usize),
    #[error("k ({k}) is larger than a counted segment (length {segment_len})")]
    KLongerThanSegment { k: usize, segment_len: usize },
}

/// Checkpoint directory I/O and integrity.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum CheckpointError {
    #[error("strict checkpoint verification failed for `{path}`: expected SHA-256 {expected}, got {actual}")]
    DigestMismatch {
        path: std::path::PathBuf,
        expected: String,
        actual: String,
    },
    #[error("checkpoint JSON error: {0}")]
    Json(#[from] serde_json::Error),
    #[error("strict checkpoint resume requires manifest at `{0}`")]
    StrictManifestMissing(std::path::PathBuf),
    #[error("graph checkpoint could not be applied: {0}")]
    InvalidGraph(String),
    #[error("checkpoint I/O: {0}")]
    Io(#[from] std::io::Error),
}
