# Recorded dependency lock and restoration procedure

Analytical baseline: R 4.5.1; TCGAbiolinks 2.38.0, edgeR 4.8.2, limma 3.66.0,
ggplot2 4.0.3, ggrepel 0.9.8 and scales 1.4.0.
jsonlite 2.0.0 supports the separate metadata-cache builder.
yaml 2.3.12 supports release metadata QA.
These values were observed in the archived/installed metadata, not invented.

renv.lock was generated from installed package descriptions with verified renv
1.3.1, including dependency closure and repository/source records.
Packaging renv was installed only in an isolated task library from the CRAN
source archive recorded in packaging_tool.json; it did not replace analysis packages.
recorded_packages.csv remains the historical snapshot with machine paths removed.
The lock is a restoration specification, not evidence that clean restoration passed.

All 103 pinned source endpoints resolved in bounded checks on 31 October 2026.
Three exact-version small packages installed into an isolated library, using existing
frozen dependencies. This is partial installation evidence, not a clean restore.
See [restoration status](RESTORATION_STATUS.md) and [source table](pinned_source_availability.csv).

Before restoration, install renv 1.3.1 into a separate tool library using the exact
source URL and SHA256 in packaging_tool.json. Then use an unused library path:

~~~sh
Rscript --vanilla tools/restore_environment.R --library ../../private_R_library_rc3 --execute
~~~

The utility requires R 4.5.1 and renv 1.3.1 and invokes restore with the lock.
Set renv state/cache paths to private writable storage. Inspect repository/source
availability before installation: old R-universe/CRAN builds may require archives.
Only scales, ggrepel and jsonlite were freshly installed into an isolated test library; historical analysis libraries were not modified.

Successful namespace/synthetic-graphics smoke tests in the existing environment
are reported separately from clean installation, Linux CI and full-data execution.
Several archived packages were built under R 4.5.2 while recorded execution used
R 4.5.1; these warnings are preserved rather than silently changing versions.
