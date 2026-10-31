# Frozen validation evidence

Amended V4 v03 is the final authority. Mandatory frozen V0-V4 scope is COMPLETE;
scientific interpretation is QUALIFIED. The canonical narrative is
[the unchanged final report](../docs/VALIDATION_REPORT.md).

V4_candidate_status_public.csv contains all 2,120 gene/model/test rows in original
order with original numeric strings. Five fields are omitted: baseline_evidence_path,
V2_worst_removed_patient_id, V3_evidence_path, prior_evidence_path and source_V2.
The remaining cells are exact archived values. The omitted original fields remain
in the source archive, not this proposed repository.

V4_claim_evidence_assessment.csv explicitly distinguishes estimable baseline
nonsignificance from unavailable evidence. The 1,272 original umbrella-category
records are not 1,272 missing checks. The final totals remain 745 supported,
103 sensitivity-qualified and 1,272 baseline-nonsignificant claims.

Completion, closure, cohort and classification-change CSVs replace machine-specific
absolute path prefixes with ARCHIVE_WORK/ or ARCHIVE_DESKTOP/ lineage aliases.
These aliases are private archive locators, not bundled downloadable files.
Their scientific values are unchanged. The other six copied validation files
are byte-identical to v03. V4_completion_summary.json is an explicitly selected
subset of the archived execution manifest.

The narrative refers to original source filenames. The raw V4_candidate_status.csv,
V4_changed_fields.csv, patient-level V2/V3 records and absolute-path manifests
stay in the private source archive. Every included transformation is recorded in
[the provenance manifest](../SOURCE_PROVENANCE_MANIFEST.csv).
Original V4 and stopped v02 are historical, not current release authorities.
No normalization, model fitting, P-value/BH calculation or candidate reselection
occurred during candidate assembly.
