# Publication instructions — future user-authorized action

This task stops at a local package. No remote repository, commit, push, upload or
publication has occurred. These steps describe a later, separately authorized action.

1. Review FINAL_RELEASE_ACCEPTANCE.md and RELEASE_BLOCKERS.md, retaining the QUALIFIED scientific interpretation and all NOT_TESTED disclosures.
2. Publish only the contents of `repository/`. Never add the enclosing candidate directory or `local_only/`. Inspect TREE.txt, LICENSE_SCOPE.csv and FILE_INVENTORY.csv before staging.
3. Verify every SHA256SUMS entry and run `Rscript --vanilla tools/run_offline_tests.R .` from the repository root. Review the license scope and TCGA/GDC attribution. Keep complete LICENSES texts.
4. Create the selected GitLab/GitHub repository only with user authorization. Retain .gitignore. Review the staged file list for private caches, RDS objects, raw counts, clinical/patient-level files, machine paths, credentials and historical notebook drafts.
5. GitLab: the default offline job is fail-closed; remote execution is NOT_TESTED here. GitHub: there is no GitHub Actions workflow. Add/review a platform workflow before claiming automated checks there. Do not include the optional real-data job by default.
6. Keep the v0.1.0-rc3 tag as a candidate. Do not label a clean restore or independent validation as passed. A later stable v0.1.0 requires a reviewed metadata/version update, refreshed manifests and final user decision.
7. Once an actual public URL/DOI exists, add it to citation metadata and refresh manifests before tagging the corresponding revision. Do not invent persistent identifiers.

Public data manifests identify GDC source files, without clinical patient mappings.
Full historical reproduction additionally requires the private clinical/export snapshot;
it is not bundled and should not be uploaded with the project. Check original data
access/acknowledgment conditions before obtaining any inputs. No new statistical
work is required or authorized by these publication instructions.
