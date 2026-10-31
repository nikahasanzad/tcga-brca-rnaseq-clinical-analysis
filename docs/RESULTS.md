# Primary and integrated results

## Separate primary analyses

Gene-wise associations were estimated using the frozen limma-voom workflow. Age remained continuous; pathological stage was a four-level factor. Statistical significance was defined as BH-FDR < 0.05 in each original model/test family.

| Analysis | Patients | Tested genes | Significant genes |
|---|---:|---:|---:|
| Primary Age | 1,072 | 25,551 | 9,063 |
| Primary Stage omnibus | 1,063 | 31,157 | 1,285 |
| Joint Age adjusted for Stage | 1,048 | 31,285 | 9,413 |
| Joint Stage omnibus adjusted for Age | 1,048 | 31,285 | 1,201 |

The Stage omnibus test detects expression heterogeneity across stage groups and has no single direction. The six pairwise contrasts characterize explicitly defined differences. Omnibus significance is not substituted for a contrast-specific FDR or stability claim.

Primary Stage groups contained 179 Stage I, 615 Stage II, 249 Stage III and 20 Stage IV tumors. The small Stage IV group warrants caution when interpreting extreme coefficients. Larger discovery counts for Age do not establish a greater biological influence of age.

Figures 1-3 illustrate the Age analysis and Figures 4-6 the Stage analysis. Selected-gene correlations, ordinary trend lines and representative expression plots are descriptive displays of the same cohort, not independent replication.

Evidence: E02-E05 and the saved Joint result tables identified in EVIDENCE_INDEX.md. Counts were read from archived results; no new P values or FDR values were computed for this package.



## Integration and the mutually adjusted model

Direct integration used 25,551 genes tested in both primary analyses. Original model-family FDR values were retained. This common display universe is distinct from the 31,157-gene Stage family and 31,285-gene Joint family.

| Association category | Primary models | Joint model |
|---|---:|---:|
| Age only | 8,550 | 7,860 |
| Stage only | 653 | 604 |
| Both | 513 | 479 |
| Neither | 15,835 | 16,608 |
| Common universe | 25,551 | 25,551 |

Of the 513 primary shared genes, 418 remained significant for both variables in the Joint model. The remaining 95 became Age only (63), Stage only (25) or Neither (7). Another 61 genes entered the Joint Both category from other primary categories, accounting for its total of 479.

Primary and Joint coefficient estimates were highly concordant genome-wide. The final Figure 8 metric table reports Pearson correlations of approximately 0.9979 for Age, 0.9948 for Stage II/I, 0.9960 for Stage III/I and 0.9994 for Stage IV/I. Use the saved full-precision table for exact values.

These comparisons also differ in cohort and preprocessing. Classification changes cannot be attributed exclusively to mutual adjustment. High global correlation is not equivalence of every gene estimate and does not establish causal independence or precise Stage IV inference.

Evidence: E06-E07. Figure 7 displays the primary categories; Figure 8 displays coefficient concordance. Historical 31,154-gene Stage comparison summaries must not be mixed with the final 25,551-gene figure universe.



## Candidate prioritization and final qualification

The historical 418 shared retained genes were ranked by the less favorable of their Joint Age and Stage FDR ranks, followed by rank sum and Ensembl ID. The top 30 were selected for a descriptive coefficient heatmap. Their identifier set matches the frozen historical-30 panel.

The historical term stable shared refers to retained significance under those primary and Joint specifications. It is not a blanket claim of independence from subtype, insensitivity to every deletion, or independent biomarker validation.

Final amended V4 v03 records completion of the frozen mandatory V0-V4 scope with scientific interpretation QUALIFIED. All 2,120 claims have their applicable frozen-scope checks. The Stage supplement filled the former 120-claim gap; the other 2,000 records were preserved.

Across the full framework, 745 claims are supported across specified checks and 103 are sensitivity-qualified. The remaining 1,272 retain the historical umbrella category but are estimable baseline-nonsignificant positive-association claims, not missing mandatory checks or computational failures.

Subtype sensitivity limits claims of independence: in the paired same-patient historical-panel Stage comparison, significance fell from 28 base-model genes to one after subtype adjustment. Individual Stage contrasts also retain influence qualifications. None of the historical 30 is baseline-significant for IV/II or IV/III.

The exploratory 12-gene panel is separate from the historical 30. Two genes overlap, giving 40 unique genes. Panel views must not be added as independent candidates. No independent cohort, mechanism, regulation, treatment benefit or clinical utility is established.

Evidence: E08-E10. Figure 9 is a prioritization display. Final claims must identify the model, test and applicable qualification; a supported omnibus does not imply a supported Stage IV contrast.
