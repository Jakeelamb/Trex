//! Paired-end **R1**/**R2** parity: same count, same ordered **pair identity** (**Phase-1 pair parity**).
//!
//! Pair identity is the FASTQ first-header-field string (stored as [`Read::id`](crate::illumina::read::Read))
//! with an optional Illumina **`/1`** or **`/2`** suffix stripped for comparison.

use crate::error::IngestError;
use crate::illumina::read::Read;

/// Strip optional **`/1`** / **`/2`** suffix on Illumina FASTQ ids so **R1** and **R2** match.
pub fn pair_identity(id: &str) -> &str {
    id.strip_suffix("/2")
        .or_else(|| id.strip_suffix("/1"))
        .unwrap_or(id)
}

pub fn validate_pair_parity(r1: &[Read], r2: &[Read]) -> Result<(), IngestError> {
    if r1.len() != r2.len() {
        return Err(IngestError::PairCountMismatch {
            r1: r1.len(),
            r2: r2.len(),
        });
    }
    for (i, (a, b)) in r1.iter().zip(r2.iter()).enumerate() {
        if pair_identity(&a.id) != pair_identity(&b.id) {
            return Err(IngestError::PairIdMismatch {
                index: i,
                r1_id: a.id.clone(),
                r2_id: b.id.clone(),
            });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn illumina_suffix_stripped() {
        assert_eq!(pair_identity("SRR/1"), "SRR");
        assert_eq!(pair_identity("SRR/2"), "SRR");
    }
}
