# Decision history and validation status

## Reading the decision history

The source notebook records an evolving analysis. Suggestions, interim checkpoints, executed methods, candidate interpretations and Results drafts occur together. The present documentation uses v09_final01 as the primary numerical baseline, v09_final02 as repeat evidence and amended V4 v03 as the current validation synthesis.

| Historical choice | Current reading |
|---|---|
| One RNA observation per participant | Avoid pseudoreplication through an explicit manifest |
| Preserve the raw count matrix | Reconstruct each analysis from the canonical input |
| Screen before excluding | Tail flags alone are not exclusion criteria |
| Separate Age and Stage filters | Match each family to its cohort and design |
| Re-estimate after six exclusions | Use preprocessing appropriate to the final cohort |
| Separate omnibus and contrasts | Overall heterogeneity has no single signed effect |
| Prioritize 30 from 418 | Descriptive post-selection, not clinical validation |

The companion DECISION_LOG.md summarizes all 62 explicitly numbered labels found in the notebook. The standalone editorial compendium preserves their source wording, the second D24 version and unnumbered documentation paragraphs. Missing labels D19, D23, D47, D53 and D65 are recorded as absent rather than invented.

An unnumbered Figure 7 decision paragraph is preserved separately. It is not relabelled as D65. Source heading repetitions are indexed without treating them as additional substantive decisions.

Source notebook: PDF pages 96-292. Decision labels are historical references, not evidence of dates, preregistration or prospective approval.



## Superseded proposals and retained limits

| Earlier notebook wording | Current correction |
|---|---|
| Numeric or ordinal Stage | Four-level factor with Stage I reference |
| Log2 described as normalization | Filtering, TMM and voom have distinct roles |
| Sex/batch-adjusted primary Age proposal | Executed primary Age-only model; sensitivities separate |
| Higher-depth or alphabetical file choice | Biospecimen and gene-assignment evidence |
| Pre-QC cohort/gene totals | Retain as historical, not final model dimensions |
| Retain missing-primary-flag exception | Frozen strict clinical inclusion rule |
| Subtype adjustment still future | Later complete-case sensitivity now exists |
| edgeR QL and Galaxy planned | No execution established in the frozen protocol |

The detailed exclusions narrative does not eliminate the current limitation on individual historical QC rationale. RIN is not a universal explanation for the six exclusions. Primary and conditional associations answer different population and adjustment questions.

Gene-function, luminal, immune, stromal and pathway explanations remain hypothesis-generating unless separately supported. Additional literature or covariate adjustment alone would not establish causality. Bulk expression cannot identify the cellular origin of a signal without additional evidence.

The historical Fisher overlap summary is retained in the compendium with its source status. It was not recomputed here and is not presented as an independent mechanistic test or new confirmatory FDR guarantee.



## Completed validation and release boundaries

| Stage | Purpose | Current recorded scope |
|---|---|---|
| V0 | Frozen protocol and panel | Reference definitions |
| V1 | Integrity and readiness | Complete with documented limitations |
| V2 | Single-patient influence | Original model paths and fixed candidates |
| V2 supplement | Three non-I primary Stage contrasts | Former 120-claim gap covered |
| V3 | Historical QC addback sensitivity | Frozen primary-model scope complete |
| Amended V4 v03 | Evidence integration | Mandatory scope COMPLETE; interpretation QUALIFIED |

The supplement covers III/II, IV/II and IV/III in the original primary Stage model. It reports 1,063 deletion cohorts, 127,560 detailed candidate records and 120 summaries. Existing Stage-versus-I and omnibus records were reconciled. These are archived execution results, not a rerun during documentation preparation.

Conditional models retain their original V3 not-applicable boundary. A completed protocol does not convert nonsignificance into a positive association or make every sensitivity flag disappear. Original V4 and the stopped amended v02 attempt remain historical provenance.

The final limitations include sparse Stage IV, incomplete individual QC rationale, weight-ID provenance, missing historical digests for supporting references and reconstructed GDC caches. Clean installation, target-runner execution, cross-platform operation and external biological validation remain unestablished.

English documentation preparation does not itself perform those checks. Repository licensing, authorship/contact information, final version/tag and distribution choices remain release decisions. No remote repository, tag or public publication was created by this task.
