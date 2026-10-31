# Data provenance and RNA sequencing selection

## Participants specimens and aliquots

The molecular source was the GDC TCGA-BRCA harmonized STAR - Counts workflow, restricted to Primary Tumor gene-expression quantification. The initial discovery set contained 1,111 files, 1,106 unique sample/vial identifiers and 1,095 participants.

Participant, sample/vial, portion, analyte, aliquot and expression file are distinct levels. More than one file for a participant is not automatically a computational duplicate. The curation aimed to retain one molecular observation per participant while preserving the identifier mapping.

| Curation stage | Participants or specimens | Files |
|---|---:|---:|
| Initial discovery | 1,095 participants | 1,111 |
| After specimen selection | 1,095 selected specimens | 1,100 |
| After aliquot selection | 1,095 participants | 1,095 |
| After strict clinical intersection | 1,093 participants | 1,093 |

Eleven participants had multiple Primary Tumor specimens. The historical selection favored the non-FFPE specimen using biospecimen annotations and available RNA integrity evidence. Alphabetical vial order was not the final selection criterion.

Five retained specimens still had paired aliquot files. Targeted expression QC of the ten files considered assigned counts, detection, STAR summary categories and pairwise expression similarity. The recorded rule favored the larger gene-assignment fraction, rather than plate identity, file size or assigned-count depth alone.

Source notebook: PDF pages 63-98. Current implementation and retained manifest remain authoritative for the exact selection. The word predefined in the notebook does not establish prospective preregistration.



## Canonical count matrix and input contract

The canonical molecular-clinical manifest identifies one selected file per retained participant and preserves file, sample and participant mappings. The count matrix was constructed from the 1,093 selected files and explicitly aligned to the clinical table.

Each parsed STAR file contained 60,664 rows and nine columns. Four summary rows were excluded from gene-level expression. The retained 60,660 Ensembl features used the unstranded raw-count field. TPM, FPKM and FPKM-UQ were not substituted for raw counts.

| Contract | Required interpretation |
|---|---|
| Raw dimensions | 60,660 gene rows by 1,093 participants |
| Observation identifiers | Unique participants with deterministic file mapping |
| Expression values | Nonnegative integer-valued raw counts |
| Gene identifiers | Consistent ordered Ensembl features |
| Sample alignment | Count columns, metadata rows and design rows match |
| Raw input object | Preserved before analysis-specific filtering |

Equality of dimensions is insufficient: identity sets and order must agree. The raw matrix remains the starting point for analysis-specific reconstruction. Completely unobserved genes were retained at the raw-matrix stage and considered during later filtering.

Assigned gene counts are not necessarily total sequenced reads. The notebook's gene-assignment fraction uses assigned counts relative to assigned and STAR-accounted summary categories; it should not be relabelled as a conventional alignment mapping rate.

Source notebook: PDF pages 79-88 and 99-120. Source descriptions of 2,676 globally zero features and aliquot-specific percentages remain historical values unless traced to their original diagnostic records.



## Input integrity and reconstruction provenance

The historical acquisition record describes exact manifest-to-disk matching and MD5 agreement for all 1,093 retained STAR files. Completeness, checksum identity and matrix alignment are distinct checks. They support input integrity, not biological validity.

The final local reproducibility lineage used reconstructed GDC caches because the original API responses were unavailable. Reconstruction used archived metadata and UUID joins. It must not be described as a newly verified live-GDC acquisition.

The retained primary outputs belong to v09_final01; v09_final02 supplies same-environment repeat evidence. A historical Stage reference object is not interchangeable with the final v09 Stage object even if its dimensions or name look similar.

This public-facing package contains documentation, selected gene-level summary tables and existing figure exports. Raw counts, clinical records, patient-to-file mappings, RData/RDS workspaces, raw logs and the full historical notebook are excluded.

GDC distinguishes open from controlled access. Expression-count access labels do not make every local derivative or clinical record suitable for unrestricted redistribution. Data access and attribution should follow the source policies. The documentation does not confer rights to redistribute third-party data.

No acquisition, count parsing, normalization or model execution was performed to prepare these documents. Source hashes and selected archived table counts were checked. Clean installation, another operating system, the intended GitLab runner and live-GDC operation remain separate verification tasks.

Source notebook: PDF pages 99-120. Current evidence: E01-E03 and E10. See DATA_AVAILABILITY.md, REPRODUCIBILITY_AND_RELEASE.md and REFERENCES.md for the distribution and evidence boundaries.
