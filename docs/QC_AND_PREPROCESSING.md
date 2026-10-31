# Quality control filtering and normalization

## Screening and sample assessment

Initial sample diagnostics included assigned gene counts, detected features, features with at least ten counts and the zero-count fraction. They characterized variation and prioritized investigation; a single distributional tail was not an automatic exclusion criterion.

Detected-feature count and zero fraction are related views of the same property. STAR fractions are also related. Counting several flags does not create several independent failures.

| Historical screening summary | Samples | Decision meaning |
|---|---:|---|
| Raw-expression distribution flags | 43 | Additional investigation |
| At least one STAR flag | 154 | Technical screening, not automatic removal |
| Raw-expression and STAR flags | 24 | Prioritized convergent candidates |
| Final historical QC exclusions | 6 | Applied before primary modelling |

Post-normalization PCA examined expression-space structure using voom log2-CPM values. PCA was exploratory QC, not the formal test of Age or Stage association. An extreme position can arise from technical or biological heterogeneity; PC1-PC2 position alone was not the exclusion rule.

The notebook provides named candidates and multivariate ranks. The current validation nevertheless retains incomplete individual exclusion reasons. RNA integrity was metadata and was not assumed to explain all six exclusions.

Source notebook: PDF pages 121-141 and 184-200. Screening totals are historical diagnostic summaries. Individual patient-level diagnostics are not included in the public package.



## Design-aware filtering TMM and voom

Low-expression filtering was performed using edgeR filterByExpr with the design for each eligible analysis. The manual CPM >= 1 in at least 20 samples rule was explored as a stringency comparison and was not adopted as the final filter.

The preliminary Age cohort retained 25,627 genes across 1,078 patients; preliminary Stage retained 31,445 across 1,069. After six exclusions, the corresponding final dimensions were 25,551 by 1,072 and 31,157 by 1,063. These processing stages must not be merged.

| Component | Output | Role |
|---|---|---|
| filterByExpr | Design-specific retained features | Defines the tested family |
| TMM | Sample normalization factors | Library-composition scaling |
| voom | log2-CPM expression and precision weights | Mean-variance-aware modelling inputs |
| limma with eBayes | Coefficients and moderated inference | Gene-wise association testing |

TMM factors describe compositional adjustments; they are not sample-quality scores and do not repair RNA degradation or other severe technical problems. Effective library size combines the analysis-specific library size and normalization factor.

Voom estimates the mean-variance relationship and observation-level precision weights. Its expression matrix supports plots and PCA; its weights support inference. A separate arbitrary log2(count + 1) transformation is not applied after voom.

Using edgeR for count handling, filtering and TMM does not mean an edgeR quasi-likelihood model was fitted. The executed primary framework was limma-voom. See REFERENCES.md for the primary TMM, voom and limma papers.



## Final reconstruction and sensitivity boundaries

After the six historical exclusions, final cohorts were reconstructed from the canonical raw matrix. Filtering, TMM and design-specific voom were re-estimated, followed by a final PCA sanity check. Preliminary weights or factors were not simply carried into the final cohort.

| Final model | Patients | Tested genes |
|---|---:|---:|
| Age | 1,072 | 25,551 |
| Stage | 1,063 | 31,157 |
| Joint | 1,048 | 31,285 |

V2 single-patient influence diagnostics hold archived expression, normalization inputs and weights fixed. V3 QC addback sensitivity changes the cohort and recalculates TMM and design-specific voom. These checks address different questions and do not substitute for one another.

Effect flags use an absolute coefficient change greater than max(0.10, 0.25 times the absolute comparator coefficient). A material sign-reversal flag also requires at least one absolute coefficient of 0.10 or more. FDR crossings are recorded separately. Lack of a flag is not an equivalence test.

Marginal 95% intervals are not simultaneous intervals. Post-selection diagnostics do not create an independent confirmatory cohort or a new conjunction-FDR guarantee. Threshold crossing alone is not proof of biological instability.

Weight matrices lack independent dimnames in the archived lineage; alignment rests on preserved-object structure, producer provenance and numerical reproduction. Historical hashes missing from the frozen reference import cannot be reconstructed retroactively by hashing a file today.

Source notebook: PDF pages 142-200 and 270-273. Current validation evidence: E10. This documentation task did not rerun any preprocessing or statistical analysis.
