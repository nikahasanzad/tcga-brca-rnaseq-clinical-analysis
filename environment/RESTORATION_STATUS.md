# Bounded dependency verification — 31 October 2026

PASS: all 103 exact pinned source URLs returned HTTP success and a gzip prefix,
including all six direct analytical dependencies. This is endpoint availability,
not validation of every archive body or proof that renv will select every endpoint.
See [source table](pinned_source_availability.csv); some Bioconductor sources resolve
through recorded exact Git commits. The lock remains byte-identical to rc2.

PASS: archive DESCRIPTION Package/Version matched scales 1.4.0, ggrepel 0.9.8,
jsonlite 2.0.0, TCGAbiolinks 2.38.0, edgeR 4.8.2 and limma 3.66.0.
ggplot2 4.0.3 archive DESCRIPTION was NOT_TESTED because the bounded 1 MB stream
limit was reached. Its pinned endpoint resolved and installed recorded namespace
version was checked; those do not substitute for archive-version verification.

PASS: scales 1.4.0, ggrepel 0.9.8 and jsonlite 2.0.0 were freshly installed from
three verified archives (1,536,246 bytes total) into a separate temporary library,
then loaded in fresh R sessions with exact versions and library location verified.
These installations used existing frozen dependencies. Historical library package
metadata was unchanged before/after; no historical package was replaced or upgraded.
[Installation evidence](../verification/bounded_installation.json) records hashes.

NOT_TESTED: a clean restoration of all 103 packages, complete transitive/system
requirements, cross-platform builds, container execution and remote CI.
The previous current-version mismatch warning is now historical; it is not evidence
that the exact old sources are unavailable. Availability checks do not prove complete
installation or full analytical reproducibility. No TCGA dataset was downloaded.
Several historical packages report builds under R 4.5.2 while execution R was 4.5.1;
those warnings are retained and no versions were changed to suppress them.
