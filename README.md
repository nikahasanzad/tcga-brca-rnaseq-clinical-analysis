# TCGA BRCA clinical association analysis

Nika Hasanzad | Version 0.1.0-rc3 | Local portfolio release candidate.

This project asks how bulk Primary Tumor gene expression is associated with age
at diagnosis and pathological stage in TCGA-BRCA. It uses one curated molecular
observation per participant, design-aware filtering, TMM, voom and limma/eBayes.
Primary Age/Stage models are distinct from the secondary mutually adjusted Joint model.

| Model | Patients | Tested genes | BH FDR below 0.05 |
|---|---:|---:|---:|
| Primary Age | 1072 | 25551 | 9063 |
| Primary Stage omnibus | 1063 | 31157 | 1285 |
| Joint Age adjusted for Stage | 1048 | 31285 | 9413 |
| Joint Stage omnibus adjusted for Age | 1048 | 31285 | 1201 |

The common 25,551-gene universe contains 513 primary Both genes; 418 retain both
associations in Joint. The historical prioritized 30 genes and exploratory 12
are separate panels with two shared genes. Final V4 comprises 745 supported,
103 sensitivity-qualified and 1,272 baseline-nonsignificant claims out of 2,120.

Frozen scientific validation is COMPLETE within bounded V0-V4 scope; interpretation
is QUALIFIED. Associations do not establish causation, regulation or clinical utility.
Stage IV contrasts are influence-sensitive; omnibus significance does not establish
pairwise significance. Subtype adjustment attenuates historical Stage findings.
Patient-specific rationales for six QC exclusions remain incomplete.
No independent biological or clinical validation has been performed.

Start with [final acceptance](FINAL_RELEASE_ACCEPTANCE.md), [Methods](docs/METHODS.md),
[Results](docs/RESULTS.md), [final validation](docs/VALIDATION_REPORT.md),
[figures](docs/FIGURE_GALLERY.html), [offline tests](tests/README.md),
[inputs](docs/INPUT_DATA.md), [execution](docs/EXECUTION.md),
[environment](environment/README.md) and [publication checklist](PUBLICATION_CHECKLIST.md).

~~~sh
Rscript --vanilla tools/run_offline_tests.R .
~~~

The offline suite is self-contained and needs no TCGA data. Full-data reproduction
additionally requires the private historical clinical/export snapshot. Public download
manifests omit patient mappings. See [data conditions](docs/IDENTIFIERS_AND_DATA_USE.md).

Primary authority is v09_final01; v09_final02 is same-environment repeat evidence;
final synthesis is amended V4 v03. Historical scripts and notebook drafts are archival.
Original analytical code, saved figures and numerical results remain unchanged.

The scoped original-code MIT and original-documentation/figure CC BY 4.0 grants
are approved; see [licensing scope](LICENSING_SCOPE.md) and [citation](CITATION.cff).
All 103 pinned source endpoints resolved. Three small exact-version packages were
installed and loaded in a separate library using existing frozen dependencies.
A full clean restoration, remote CI and full-data reproduction remain NOT_TESTED.
See [restoration evidence](environment/RESTORATION_STATUS.md),
[remaining requirements](RELEASE_BLOCKERS.md) and [publication instructions](PUBLICATION_INSTRUCTIONS.md).
This is READY FOR USER-AUTHORIZED PUBLICATION as a qualified research portfolio,
with the disclosed limits. No repository was created remotely or published.
