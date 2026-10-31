# Figure guide and interpretation

## Age figures

The nine copied PNG exports come from v09_final01 and are unchanged. Each has a stable repository filename, an English legend and a source checksum. The local gallery displays the actual figures with the corresponding interpretation notes.

| Figure | Display | Interpretation |
|---|---|---|
| 1 | Age coefficient per decade versus BH-FDR | Continuous association, not young-versus-old DE |
| 2 | Selected gene expression versus age | Weighted inference with descriptive Pearson/OLS overlays |
| 3 | Largest significant protein-coding Age effects | Extreme-effect subset, not the shared historical 30 |

Figure 2 illustrates ESR1, SOX10, GABRP and KRT5. Primary coefficients and FDR values come from the weighted limma-voom model. Pearson correlations and ordinary trend lines are descriptive summaries on selected genes. Free y-axis scales support within-gene reading; visual slope angles do not compare effects across panels.

Strong statistical evidence can coexist with modest individual-level correlations and substantial between-tumor heterogeneity. These plots do not demonstrate clinical predictive performance.

Figure 3 uses an effect-size ranking among significant protein-coding genes. It is a different selection from the 30 shared genes in Figure 9. Directional imbalance in an extreme-effect subset is not a transcriptome-wide directional conclusion.

Source notebook: PDF pages 217-221. Read the full figure legends with the existing export. No figure was regenerated, relabelled or statistically re-estimated for this documentation package.



## Stage figures

| Figure | Display | Scale and hypothesis |
|---|---|---|
| 4 | Transcriptome-wide Stage overview | Unsigned maximum coefficient; omnibus FDR |
| 5 | Selected expression distributions across stages | Descriptive boxplots and unlabeled tumor points |
| 6 | Stage-versus-I coefficient heatmap | Raw coefficients, symmetric color cap, no row z scores |

Figure 4 separates statistical evidence for overall Stage association from the magnitude of an explicit Stage-versus-I coefficient. The omnibus has no single positive or negative direction.

Figure 5 shows ACKR1, ADH1B, ADIPOQ and PIR. Group sizes are Stage I 179, II 615, III 249 and IV 20. Its FDR annotations refer to the primary omnibus test. Unequal group sizes and free panel scales must remain explicit.

Figure 6 displays the strongest omnibus-associated protein-coding genes. Raw coefficients retain direction and magnitude relative to Stage I. Colors saturate at the stated cap of +/-2.5; the underlying coefficients are not truncated. This encoding differs from the column scaling in Figure 9.

Nonmonotonic patterns support a descriptive categorical representation but do not prove the absence of every possible ordered component. Stage groups are cross-sectional categories, not a within-patient progression series. Significant-gene counts across contrasts are not biological distances.

Source notebook: PDF pages 241-250. Stage-IV-specific effect interpretations must retain the final influence qualifications, even when an omnibus result or whole-genome correlation looks stable.



## Integrated figures and final checks

| Figure | Display | Required qualification |
|---|---|---|
| 7 | Four primary association categories | Denominator 25,551; statistical categories |
| 8 | Primary versus Joint coefficients | Cohort and preprocessing also differ |
| 9 | Top 30 shared retained candidates | Selection-specific stability; independent column scales |

Figure 7 preserves the common tested denominator. Age only and Stage only describe thresholded statistical membership, not biological specificity.

Figure 8 uses 25,551 genes per panel. The saved final annotations round Stage II/I to 0.9948 and Stage III/I to 0.9960; some notebook drafts use slightly different rounding. High global concordance does not establish equivalence of individual estimates or precision for Stage IV.

Figure 9 uses 30 genes selected from 418 primary-to-Joint retained shared associations. The Stage selection criterion is omnibus significance, not significance in every individual Stage contrast. Each column is divided by its own maximum absolute coefficient among the 30 genes; colors therefore cannot compare absolute magnitudes across columns.

The phrase stable shared in the unchanged historical figure title is limited to the selection definition. It is not an independent validation claim. Final subtype and influence qualifications remain essential.

The gallery and FIGURE_LEGENDS.md give all nine legends. The editorial compendium additionally preserves historical caption variants, the unnumbered Figure 7 decision and the revised Figure 9 wording on pages 291-292. Visual review of the nine PNGs found readable English labels at review size; final journal-specific sizing remains a separate delivery check.
