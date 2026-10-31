# Final release acceptance — 0.1.0-rc3

**READY FOR USER-AUTHORIZED PUBLICATION** as a qualified research portfolio.
Date: 31 October 2026. Author/copyright holder for authorized original material:
Nika Hasanzad. This is a local release candidate; no online publication occurred.

| Status | Acceptance item | Evidence and boundary |
|---|---|---|
| PASS | Approved scoped licensing | Complete MIT and CC BY 4.0 legal texts; LICENSING_SCOPE.md and LICENSE_SCOPE.csv; no grant over TCGA or third-party resources. Ownership relies on author approval and audited lineage. |
| PASS | Citation | CFF 1.2.0 schema, verified author name, MIT software license; original docs/figures CC BY separately scoped; no invented DOI/contact. |
| PASS | Source integrity | rc2 source manifests reverified; rc1, rc2 and historical files hash-preserved; current per-file inventory and SHA256SUMS. |
| PASS | Frozen science | Original R pipeline, results, selected figures, V4 evidence, scientific fixture and unchanged renv.lock retained by hash; no new analysis. |
| PASS | Offline QA | Static/structural, synthetic runtime, optional-QC, public-table, help, integrity and synthetic plotting checks rerun; namespace smoke, local links, fixtures, figures, safety and ignore checks pass. |
| PASS | Dependency endpoints | All 103 pinned source endpoints responded with gzip headers; this is not full archive validation. |
| PASS | Bounded installations | scales 1.4.0, ggrepel 0.9.8, jsonlite 2.0.0 installed/loaded in a new library using existing dependencies; historical libraries unchanged. |
| PASS | Distribution boundary | Only repository/ is public; no raw data, clinical tables, RDS objects, notebook drafts, credentials or machine paths in public tree. Fixed scoped TCGA provenance labels retained with data-use notice. |
| PASS | Local CI structure | GitLab YAML, fail-closed offline job and pinned image digest validated; optional real-data job remains separate/manual. |
| NOT_TESTED | Complete clean restoration and system dependency closure | 103-package fresh restore not executed; source endpoint availability and three installs do not establish it. |
| NOT_TESTED | Remote CI/container execution | No remote Runner or Docker run occurred. |
| NOT_TESTED | Full-data / independent validation | No full-data run; private historical input snapshots are not bundled. No independent biological or clinical validation. |
| NEEDS_REVISION | Optional platform/scientific documentation follow-up | GitHub Actions absent if that platform is selected; six historical per-patient QC rationales remain incomplete. See RELEASE_BLOCKERS.md. |

The licensing/source-availability release blockers from rc2 are closed within the
explicit verification bounds. The remaining items do not block a transparently
qualified portfolio release; they must block any stronger claim of clean installation,
independent reproduction, complete QC rationale or validated clinical performance.
Scientific validation remains COMPLETE only within bounded V0–V4 and interpretation
remains QUALIFIED. Frozen model results are not altered by this release work.

Consult RELEASE_QA_MATRIX.csv, verification/, environment/RESTORATION_STATUS.md,
RELEASE_BLOCKERS.md and PUBLICATION_INSTRUCTIONS.md. Publish only repository/ after
separate authorization. Stop point: local package and acceptance report.
