# Distribution boundaries

repository/ is the proposed public tree. local_only/ is never published.
Raw count files, raw/curated clinical or biospecimen data, patient/file mappings,
individual deletion diagnostics, RDS/RData caches and workspaces, raw logs,
historical notebook drafts and machine-specific provenance remain local.

Public file manifests contain archived open-file download metadata, not patient
phenotypes. Frozen TCGA labels in code and the text syntax fixture are assessed
in IDENTIFIERS_AND_DATA_USE.md. The fixture contains language expressions only.
It is not a patient-data object and eliminates the private RDS test dependency.

.gitignore is a convenience, not an access-control boundary. Inspect a future
staged file list using PUBLICATION_CHECKLIST.md before any commit or upload.
Do not publish the enclosing candidate folder or an archive containing local_only.
Original project licensing is approved within LICENSING_SCOPE.md.
