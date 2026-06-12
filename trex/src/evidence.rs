//! Typed assembly evidence records shared by graph stages and benchmark artifacts.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

pub const EVIDENCE_LEDGER_SCHEMA_VERSION: u64 = 1;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EvidenceKind {
    MateBridgeExistingEdge,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EvidenceSourceStage {
    Phase2MateBridge,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct SupportCounts {
    pub observed: u64,
    pub eligible: u64,
    pub supporting: u64,
    pub applied: u64,
}

impl SupportCounts {
    pub fn confidence(&self) -> ConfidenceScore {
        ConfidenceScore::from_ratio(self.supporting, self.eligible)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ConfidenceScore {
    pub numerator: u64,
    pub denominator: u64,
    pub value: f64,
}

impl ConfidenceScore {
    pub fn from_ratio(numerator: u64, denominator: u64) -> Self {
        let value = if denominator == 0 {
            0.0
        } else {
            numerator as f64 / denominator as f64
        };
        Self {
            numerator,
            denominator,
            value,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EvidenceRecord {
    pub kind: EvidenceKind,
    pub source_stage: EvidenceSourceStage,
    pub support: SupportCounts,
    pub confidence: ConfidenceScore,
    pub counters: BTreeMap<String, u64>,
}

impl EvidenceRecord {
    pub fn new(
        kind: EvidenceKind,
        source_stage: EvidenceSourceStage,
        support: SupportCounts,
    ) -> Self {
        Self {
            kind,
            source_stage,
            support,
            confidence: support.confidence(),
            counters: BTreeMap::new(),
        }
    }

    pub fn with_counter(mut self, name: impl Into<String>, value: u64) -> Self {
        self.counters.insert(name.into(), value);
        self
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EvidenceLedger {
    pub schema_version: u64,
    pub records: Vec<EvidenceRecord>,
}

impl EvidenceLedger {
    pub fn new() -> Self {
        Self {
            schema_version: EVIDENCE_LEDGER_SCHEMA_VERSION,
            records: Vec::new(),
        }
    }

    pub fn from_records(records: Vec<EvidenceRecord>) -> Self {
        Self {
            schema_version: EVIDENCE_LEDGER_SCHEMA_VERSION,
            records,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub fn push(&mut self, record: EvidenceRecord) {
        self.records.push(record);
    }
}

impl Default for EvidenceLedger {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::{
        EvidenceKind, EvidenceLedger, EvidenceRecord, EvidenceSourceStage, SupportCounts,
        EVIDENCE_LEDGER_SCHEMA_VERSION,
    };

    #[test]
    fn evidence_record_construction_computes_confidence() {
        let record = EvidenceRecord::new(
            EvidenceKind::MateBridgeExistingEdge,
            EvidenceSourceStage::Phase2MateBridge,
            SupportCounts {
                observed: 10,
                eligible: 7,
                supporting: 4,
                applied: 3,
            },
        )
        .with_counter("trusted_endpoint_pairs", 5);

        assert_eq!(record.confidence.numerator, 4);
        assert_eq!(record.confidence.denominator, 7);
        assert!((record.confidence.value - (4.0 / 7.0)).abs() < f64::EPSILON);
        assert_eq!(record.counters["trusted_endpoint_pairs"], 5);
    }

    #[test]
    fn evidence_ledger_defaults_to_current_schema() {
        let ledger = EvidenceLedger::default();
        assert_eq!(ledger.schema_version, EVIDENCE_LEDGER_SCHEMA_VERSION);
        assert!(ledger.is_empty());
    }
}
