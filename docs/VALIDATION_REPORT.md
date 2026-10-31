# Final Validation Evidence Integration — TCGA-BRCA, amended V4 v03

Date: 29 October 2026. **Frozen V0–V4 protocol: COMPLETE within its bounded mandatory scope. Final scientific decision: QUALIFIED.**

## Executive summary

The verified V2 Stage Contrast Coverage Supplement closes the only previously recorded mandatory coverage gap: `III_vs_II`, `IV_vs_II`, and `IV_vs_III` in the original primary Stage model. All 120 corresponding claim records now contain V2 influence evidence. The other 2,000 records retain their original values and raw CSV record bytes, including classification, evidence, numeric precision and quoting. The full framework remains 2,120 unique gene/model/test records in the original order.

This amendment integrates archived results only. It performs no normalization, voom estimation, statistical fitting, P-value generation, BH calculation, filtering, candidate selection, edgeR quasi-likelihood, external validation or publication. All historical V0–V4 evidence and the supplement remain unchanged. The frozen historical 30-gene panel and exploratory 12-gene panel remain separate, with two shared genes and a 40-gene union.

Across the complete framework: **745 supported**, **103 sensitivity-qualified**, and **1272 in the original insufficient-information/not-estimable positive-claim category**. Every record in that last category is presently estimable and baseline-nonsignificant; none represents a remaining missing mandatory check or computationally non-estimable baseline. The companion `V4_claim_evidence_assessment.csv` makes this distinction explicit without altering the preserved records.

Completion describes execution and coverage of the frozen protocol. It does not assert that every association is significant, influence-free, independent of subtype or independently validated. Sensitivity and nonsignificance are valid outcomes rather than reasons to select replacement candidates or expand the protocol.

## Verification and bounded amendment

The supplement final report, all 39 indexed supplement outputs, original V4 outputs and original V2/V3 final manifests were verified against archived SHA256 and size. The initial v02 preflight stopped before integration because one Finder `.DS_Store` file differed from historical snapshots. That stopped attempt is preserved. This v03 amendment explicitly treats that single OS directory-layout metadata file as non-scientific: both mismatching historical references are recorded in `V4_non_scientific_metadata_exceptions.csv`. No scientific data, code, protocol, baseline, report or candidate artifact is exempted. The exception does not establish who or what changed the Finder metadata. The supplement's 120 unique summaries, all archived mandatory/numerical/coverage checks, 1,063 deletion reconciliation records, environment agreement and empty non-estimable-case table were checked before integration. No archived statistical script was executed.

The supplement used fixed primary Stage expression, normalization inputs and voom weights for 1,063 deletion cohorts, with 31,157 genes per contrast, before candidate extraction. It reported 127,560 detailed candidate/deletion records and 120 summaries. Full-cohort contrast checks had maximum absolute difference 2.22045e-16; four prespecified direct-refit checks had maximum difference 6.63913e-14. Existing Stage-versus-I and omnibus results reconciled over 170,080 archived records. These are verified prior execution records, not fits rerun in this amendment.

The authorized changes are restricted to V2 coverage/count/flag/range/influential-patient/source fields and their resulting classification/reason codes for the 120 keys. Their baseline, V3, prior sensitivity, panel, estimand and claim-limit fields are preserved. `V4_changed_fields.csv` records every changed cell. `V4_classification_changes.csv` includes all 120 claims even where the category remains unchanged; `V4_unchanged_2000_verification.csv` proves preservation of the remaining records.

| Old classification | New classification | Claims |
|---|---|---:|
| insufficient information / not estimable | supported across the specified checks | 17 |
| insufficient information / not estimable | insufficient information / not estimable | 99 |
| insufficient information / not estimable | sensitivity-qualified | 4 |

## Unchanged classification and statistical rules

1. Missing required evidence, unavailable baseline significance or baseline FDR >=0.05 retain the original umbrella classification **insufficient information / not estimable**. Reason codes and explicit evidence-status fields distinguish actual missingness/non-estimability from an estimable nonsignificant result. All 99 baseline-nonsignificant supplemented claims remain in this category, now with complete V2 coverage.
2. A baseline-significant claim with complete applicable evidence is **sensitivity-qualified** if any applicable V2 effect/sign flag or FDR crossing, V3 effect/sign flag or FDR crossing, or same-cohort/same-formula weighting qualification is recorded. Otherwise it is **supported across the specified checks**. Existing V3 and prior evidence continue to contribute; the supplement does not replace them.
3. Primary models require V2 plus all seven V3 cases. Restricted/conditional models require their existing V2 and paired weighting evidence; V3 remains explicitly outside their frozen scope. The completion statement does not extend their coverage.
4. The V0 definitions remain FDR <0.05 within each original full model/test gene family; marginal moderated-t 95% CIs; effect change >max(0.10,0.25 x absolute comparator effect); sign reversal with at least one absolute effect >=0.10. Age is per decade. No FDR was recalculated in this amendment.
5. Stage omnibus is a joint F test with no coefficient direction. Individual contrasts use their own moderated-t estimates and full-family FDR. Omnibus significance is never substituted for contrast significance or contrast-specific stability. Primary between-tumor and subtype-adjusted conditional associations remain different estimands.

## Supplementary Stage findings

Counts are candidate claims with at least one flagged deletion, not counts of patients; effect, sign and FDR columns can overlap. Historical/exploratory views also overlap and must not be summed.

| Contrast | Panel | Baseline significant | V2 effect flags | V2 sign flags | V2 FDR crossings |
|---|---|---:|---:|---:|---:|
| III_vs_II | union40 | 18/40 | 0 | 0 | 3 |
| III_vs_II | historical30 | 17/30 | 0 | 0 | 3 |
| III_vs_II | exploratory12 | 2/12 | 0 | 0 | 0 |
| IV_vs_II | union40 | 1/40 | 29 | 5 | 3 |
| IV_vs_II | historical30 | 0/30 | 27 | 5 | 0 |
| IV_vs_II | exploratory12 | 1/12 | 4 | 0 | 3 |
| IV_vs_III | union40 | 2/40 | 27 | 3 | 3 |
| IV_vs_III | historical30 | 0/30 | 24 | 2 | 0 |
| IV_vs_III | exploratory12 | 2/12 | 5 | 1 | 3 |

`III_vs_II` has 18/40 baseline-significant union claims (17/30 historical). It has no predeclared effect or material-sign flags, but three claims have FDR crossings. Those crossings and the existing V3 addback evidence determine individual qualifications; absence of a flag is not an equivalence test.

The Stage IV comparisons show substantial coefficient influence: 29/40 effect-flagged claims for `IV_vs_II` and 27/40 for `IV_vs_III`, with 5 and 3 material sign reversals, respectively. These diagnostics apply even to baseline-nonsignificant estimates. Stage IV has only 20 patients in the original primary Stage cohort. The largest recorded absolute coefficient movement is about 0.486 log2 units. These are leave-one-out influence observations, not longitudinal stage progression, causal effects or a justification to remove an influential patient.

None of the historical 30 is baseline-significant in either new Stage IV contrast. In the union there are only one and two baseline-significant claims, respectively; they remain exploratory-panel observations and their final status must be read with the recorded influence/addback qualifications. The historical 30's supported Stage omnibus claims do not establish robust differences between Stage IV and II/III. Existing Stage IV/I qualifications remain unchanged in the other 2,000 rows.

## Historical 30: primary and joint claim status

The final column below denotes baseline-nonsignificant positive-association claims with complete applicable evidence, not computational failures.

| Model | Test | Baseline significant | Supported | Sensitivity-qualified | Baseline nonsignificant |
|---|---|---:|---:|---:|---:|
| age | age_10y | 30 | 30 | 0 | 0 |
| stage | Stage_omnibus | 30 | 30 | 0 | 0 |
| stage | II_vs_I | 21 | 17 | 4 | 9 |
| stage | III_vs_I | 1 | 0 | 1 | 29 |
| stage | IV_vs_I | 0 | 0 | 0 | 30 |
| stage | III_vs_II | 17 | 16 | 1 | 13 |
| stage | IV_vs_II | 0 | 0 | 0 | 30 |
| stage | IV_vs_III | 0 | 0 | 0 | 30 |
| joint | age_10y | 30 | 30 | 0 | 0 |
| joint | Stage_omnibus | 30 | 30 | 0 | 0 |
| joint | II_vs_I | 23 | 21 | 2 | 7 |
| joint | III_vs_I | 2 | 1 | 1 | 28 |
| joint | IV_vs_I | 0 | 0 | 0 | 30 |

## Prior subtype, weighting and target population

The unchanged prior evidence is retained in the byte-identical `V4_prior_evidence_long.csv`, `V4_prior_adjustment_comparison.csv`, `V4_prior_weighting_comparison.csv` and `V4_cohort_coverage.csv` copies. Primary Stage and Joint Stage have different cohorts and testing families; integrated S9 common-universe FDR is not a substitute for full-model FDR.

The RIN and purity complete-case blocks contain 1,047 and 974 patients; subtype and combined blocks contain 938 and 905 women. Relative to the 1,036-woman reference cohort, subtype/combined blocks omit 98/131 women. Common-weight and design-specific paths inherit archived 1,048-person TMM inputs; they are paired sensitivity comparisons, not independent replication. ABSOLUTE purity is a matched DNA-portion measure, not a direct measurement of RNA cell composition.

In same-patient subtype models, historical-panel Stage omnibus significance fell from 28 base-model genes to one after adjustment. This limits claims of independence from subtype without automatically invalidating the primary estimand. CCL21's archived subtype-adjusted omnibus FDR is approximately 0.01875 under common weights and 0.03751 under design-specific weights; the latter crosses 0.05 in seven original V2 deletions. This remains a conditional, influence/weighting-qualified association rather than independent biomarker validation. No prior conditional classification was altered here.

The exploratory 12 remain MAGEC2, PIR, CCL21, OR56A3, SPPL2C, MAGEA6, MAGEA3, LINC02055, ZDHHC7, CD1C, NHLH1 and AC124045.1. They were selected after prior sensitivity results and must retain the exploratory label.

## Scientific and provenance limitations

- Individual reasons for the six historical QC exclusions are incomplete. The documented shared reason is concordant severe technical QC and multivariate expression-space abnormality; RIN must not be presented as the reason for all six. V3 addback is sensitivity analysis, not correction of the primary analysis.
- V3 produced no predeclared candidate effect/sign flags but did produce FDR crossings; those existing qualifications remain included. V2 fixes normalization/expression/weights for deletion influence, whereas V3 recalculates TMM and design-specific voom for cohort changes. The checks answer different questions.
- Stage IV is sparse (20 primary/joint patients, 17 subtype, 16 combined); omnibus stability cannot hide individual-contrast influence or nonsignificance.
- Marginal CIs are not simultaneous intervals. Post-selection checks are not a new confirmatory family or established conjunction-FDR guarantee; FDR crossings alone do not prove biological instability, and lack of crossings does not prove equivalence.
- Weight matrices lack independent dimnames; alignment relies on preserved-object structure, producer provenance and numerical reproduction. Twenty-five supporting references lacked historical digests imported into V0. Current preservation cannot create those missing historical checksums; the 28 nonmandatory V1 NOT_CHECKED entries remain documented.
- Original GDC API responses were unavailable and caches were reconstructed using archived metadata/UUID joins. Clean installation, intended-runner operation, live-GDC and cross-platform reproducibility have not been established by this amendment.
- No new batch/TSS adjustment, nonlinear Age model, interaction, prediction model, edgeR QL, independent cohort or external/clinical validation was added. No regulation, causality, mechanism, within-patient progression, treatment effect or clinical utility is established.

## Completion and release assessment

**Full frozen V0–V4 mandatory scope: COMPLETE. Required claim coverage: COMPLETE. Remaining mandatory coverage gaps: ZERO. Overall scientific interpretation: QUALIFIED.** All 2,120 claims have their applicable frozen-scope checks. The 120 previously uncovered primary Stage contrasts now have verified V2 plus existing V3 evidence. Conditional models' V3 N/A status is unchanged by design. Scientific qualifications and baseline nonsignificance remain in the final classifications; they do not require further fitting to declare this bounded protocol complete.

The amendment changes only the integration record and completion assessment in a new directory. It does not amend V0, overwrite original V2/V3/V4 outputs, change the panels or replace primary results. An old V4 report that records the former gap remains historically accurate for its version; this v03 report supersedes its completion assessment for release use.

Prioritized GitLab/Portfolio release checklist:

1. Approve model/test-specific release wording: preserve historical versus exploratory labels, subtype/complete-case limits, sparse Stage IV influence, baseline nonsignificance and QC provenance. Use the amended report as current status.
2. Assemble a new versioned English release package containing audited v09 code, amended report/status/change tables, evidence index and SHA256 manifests. Preserve the source archives.
3. Provide portable path/data/cache instructions and an installable environment specification. Verify installation and the intended GitLab runner separately; archived package snapshots and existing base-R/help CI checks are not proof of a fresh full-data or cross-platform run.
4. Review final figures/captions and portfolio narrative at delivery size, then select version/tag, license and changelog. Publish only when separately authorized. No publication or new verification run was performed here.

## Reproducibility and delivery

`V4_integration_checks.csv`, `V4_archived_manifest_verification.csv`, `V4_source_preservation.csv`, `V4_evidence_input_manifest.csv`, `V4_output_manifest.csv` and `V4_execution_manifest.json` document archived comparisons, exact join scope, source preservation, current output hashes and execution context. `V4_protocol_completion.csv`, `V4_coverage_closure.csv` and the empty (header-only) `V4_coverage_gaps.csv` distinguish closed coverage from scientific qualifications. `V4_claim_evidence_assessment.csv` distinguishes estimable nonsignificance from missing/non-estimable evidence.

Source-preservation results and final file counts are recorded in the execution manifest after integration. The output manifest excludes itself to avoid a self-referential hash. The integration uses Python standard-library file/hash/CSV operations only; no statistical package or archived analysis script was run.

**Stop point: amended V4 complete. No new analyses or GitLab/Portfolio publication initiated.**
