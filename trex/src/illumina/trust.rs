//! Report-only read-trust diagnostics for the k-mer count table.

#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct TrustDiagnosticsReport {
    pub schema_version: u64,
    pub report_only: bool,
    pub k: usize,
    pub trusted_threshold: u64,
    pub total_unique_kmers: usize,
    pub trusted_kmers: usize,
    pub rejected_kmers: usize,
    pub singleton_kmers: usize,
    pub trusted_fraction: f64,
    pub max_multiplicity: u64,
    pub multiplicity_n50: u64,
    pub buckets: Vec<TrustMultiplicityBucket>,
    pub correction_candidates: Vec<CorrectionCandidateSummary>,
}

#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct TrustMultiplicityBucket {
    pub label: String,
    pub min_count: u64,
    pub max_count: Option<u64>,
    pub unique_kmers: usize,
    pub trusted_kmers: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct CorrectionCandidateSummary {
    pub read_id: String,
    pub start: usize,
    pub end: usize,
    pub reason: String,
}

pub fn build_trust_diagnostics(
    k: usize,
    trusted_threshold: u64,
    counts: &[(Vec<u8>, u64)],
) -> TrustDiagnosticsReport {
    let total_unique_kmers = counts.len();
    let trusted_kmers = counts
        .iter()
        .filter(|(_, count)| *count >= trusted_threshold)
        .count();
    let singleton_kmers = counts.iter().filter(|(_, count)| *count == 1).count();
    let max_multiplicity = counts.iter().map(|(_, count)| *count).max().unwrap_or(0);
    let multiplicity_n50 = multiplicity_n50(counts);
    let rejected_kmers = total_unique_kmers.saturating_sub(trusted_kmers);
    TrustDiagnosticsReport {
        schema_version: 1,
        report_only: true,
        k,
        trusted_threshold,
        total_unique_kmers,
        trusted_kmers,
        rejected_kmers,
        singleton_kmers,
        trusted_fraction: fraction(trusted_kmers, total_unique_kmers),
        max_multiplicity,
        multiplicity_n50,
        buckets: build_buckets(trusted_threshold, counts),
        correction_candidates: Vec::new(),
    }
}

fn build_buckets(
    trusted_threshold: u64,
    counts: &[(Vec<u8>, u64)],
) -> Vec<TrustMultiplicityBucket> {
    let mut buckets = vec![
        TrustMultiplicityBucket::new("singleton", 1, Some(1)),
        TrustMultiplicityBucket::new("low_support", 2, trusted_threshold.checked_sub(1)),
        TrustMultiplicityBucket::new("trusted_low", trusted_threshold, Some(trusted_threshold)),
        TrustMultiplicityBucket::new(
            "trusted_mid",
            trusted_threshold.saturating_add(1),
            trusted_threshold.checked_mul(4),
        ),
        TrustMultiplicityBucket::new(
            "trusted_high",
            trusted_threshold.saturating_mul(4).saturating_add(1),
            None,
        ),
    ];

    for (_, count) in counts {
        for bucket in &mut buckets {
            if bucket.contains(*count) {
                bucket.unique_kmers += 1;
                if *count >= trusted_threshold {
                    bucket.trusted_kmers += 1;
                }
                break;
            }
        }
    }
    buckets
        .into_iter()
        .filter(|bucket| match bucket.max_count {
            Some(max) => bucket.min_count <= max,
            None => true,
        })
        .collect()
}

impl TrustMultiplicityBucket {
    fn new(label: &str, min_count: u64, max_count: Option<u64>) -> Self {
        Self {
            label: label.to_string(),
            min_count,
            max_count,
            unique_kmers: 0,
            trusted_kmers: 0,
        }
    }

    fn contains(&self, count: u64) -> bool {
        count >= self.min_count
            && match self.max_count {
                Some(max) => count <= max,
                None => true,
            }
    }
}

fn multiplicity_n50(counts: &[(Vec<u8>, u64)]) -> u64 {
    let total: u128 = counts.iter().map(|(_, count)| u128::from(*count)).sum();
    if total == 0 {
        return 0;
    }
    let midpoint = total.div_ceil(2);
    let mut multiplicities: Vec<u64> = counts.iter().map(|(_, count)| *count).collect();
    multiplicities.sort_unstable_by(|a, b| b.cmp(a));
    let mut running = 0u128;
    for count in multiplicities {
        running += u128::from(count);
        if running >= midpoint {
            return count;
        }
    }
    0
}

fn fraction(numerator: usize, denominator: usize) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

#[cfg(test)]
mod tests {
    use super::build_trust_diagnostics;

    #[test]
    fn trust_report_counts_threshold_strata_without_correction() {
        let counts = vec![
            (b"AAA".to_vec(), 1),
            (b"AAC".to_vec(), 2),
            (b"ACC".to_vec(), 3),
            (b"CCC".to_vec(), 12),
        ];

        let report = build_trust_diagnostics(3, 3, &counts);

        assert_eq!(report.schema_version, 1);
        assert!(report.report_only);
        assert_eq!(report.total_unique_kmers, 4);
        assert_eq!(report.trusted_kmers, 2);
        assert_eq!(report.rejected_kmers, 2);
        assert_eq!(report.singleton_kmers, 1);
        assert_eq!(report.max_multiplicity, 12);
        assert_eq!(report.multiplicity_n50, 12);
        assert!(report.correction_candidates.is_empty());
        assert!(report
            .buckets
            .iter()
            .any(|bucket| bucket.label == "low_support" && bucket.unique_kmers == 1));
    }
}
