# Data availability and distribution boundary

The original data source is the NCI Genomic Data Commons TCGA-BRCA project. Access classifications must be checked at the source; open expression-count records do not automatically authorize redistribution of every local clinical or derived object. See the [GDC data-access guide](https://gdc.cancer.gov/access-data/data-access-processes-and-tools).

This release candidate includes unchanged audited v09 code, existing English figures, selected gene-level and aggregate tables, and privacy-filtered validation evidence. It excludes raw counts, clinical tables, participant/aliquot mappings, patient-level deletion records, RData/RDS objects, raw logs, credentials and the full historical notebook.

The original GDC API responses were unavailable in the reproducibility archive. Caches were reconstructed from archived metadata and UUID joins. Recorded local reproduction is not a current live-GDC acquisition test. Input retrieval, reference snapshot, environment installation and the intended runner must be documented before distributing a fully executable release.

The private source archive retains the complete lineage. Public documentation uses project-relative evidence identifiers and fingerprints rather than personal filesystem paths. Gene identifiers and aggregate statistics do not identify individual participants, while copied scatterplots contain unlabeled cohort observations.

Source-data attribution and any data-use conditions remain those of the data providers. The author is Nika Hasanzad. Original project licensing is approved within LICENSING_SCOPE.md. Fixed code labels and archived open-file metadata are assessed in IDENTIFIERS_AND_DATA_USE.md; raw clinical data and patient-level diagnostics remain local. No rights grant, publisher acceptance, new patient-data sharing agreement or clinical validation is implied by this documentation package.
