//! One preprocessed read record (header id + **ACGTN** upper sequence).

#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct Read {
    pub id: String,
    pub sequence: Vec<u8>,
}

impl Read {
    pub fn shortest_length(reads: &[Read]) -> Option<usize> {
        reads.iter().map(|r| r.sequence.len()).min()
    }
}
