Teaching/decision artifacts for the v0.7.0 program.

- delta_null_distributions.html — "The distribution was the object": null
  delta-hat ridgelines by q, matched-null (steelblue) highlighting,
  percentile-on-matched-null reading. Data embedded from nulls.json
  (Stage 3/4 sweep, PRE-production-pipeline; regenerate both from the
  Stage 6 production-faithful draws when they land — the data block is
  the __DATA__ splice; extraction script pattern in session log
  2026-07-04).
- delta_envelope_explorer.html — earlier cut-line/envelope explorer
  (superseded by the distribution-native reading; kept for the record).
- nulls.json — per-(k,N,q) null summaries: 56-bin histogram, 50-point
  quantile grid, median/p95/p99 + bootstrap CI (5,000 draws/cell).
