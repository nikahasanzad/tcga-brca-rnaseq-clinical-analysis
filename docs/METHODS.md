# Methods

## Study and input definition

TCGA-BRCA bulk RNA-seq STAR - Counts Primary Tumor files were matched to a strictly curated patient-level clinical table. Age at diagnosis was converted from days to years. AJCC substages were collapsed to major Stages I-IV; Stage X and missing values were unavailable for stage-specific inference. The canonical raw matrix contained 60,660 features and 1,093 participants.

One molecular observation per participant was retained. Biospecimen evidence supported non-FFPE specimen selection where paired alternatives existed. Paired aliquot files were assessed using targeted expression QC and gene-assignment fraction. These are distinct curation levels. See [sample selection](DATA_PROVENANCE_AND_SELECTION.md).

## QC and analysis-specific preprocessing

Raw-count and STAR summaries supported screening. Final historical exclusions involved six samples with concordant technical and multivariate expression-space abnormalities. Individual exclusion provenance remains incomplete; the final validation retains that limitation. Filtering, TMM normalization and design-specific voom were re-estimated after exclusions.

Low-expression features were filtered using filterByExpr for the eligible cohort and design of each model. TMM and voom provide distinct composition-scaling and mean-variance/precision-weight components. The tested families were 25,551 genes for Age, 31,157 for Stage and 31,285 for Joint. See [QC and preprocessing](QC_AND_PREPROCESSING.md).

## Association models and testing families

The primary Age model was expression ~ age_years. Its coefficient was reported additionally per decade. The primary Stage model was expression ~ stage_group with Stage I as the treatment reference. Overall Stage association used the moderated omnibus F test of the three non-reference coefficients. All six pairwise contrasts were evaluated separately, with BH correction over the full retained Stage family for each contrast.

The secondary Joint model was expression ~ age_10y + stage_group in 1,048 complete-case participants. Age used a moderated coefficient test; Stage used the joint test of its three coefficients. BH-FDR < 0.05 was interpreted within each original model/test family. Omnibus significance was never a substitute for a contrast-specific FDR.

## Integration and visualization

Direct integration used the 25,551 genes tested in both primary analyses, retaining original-family FDR values. Of 513 primary Both genes, 418 remained Both in Joint. The historical top 30 were ranked by the worse of the two Joint FDR ranks, then rank sum and Ensembl ID. Figure 9 coefficients were scaled independently within columns by their maximum absolute value among those 30 genes.

Primary-to-Joint comparisons differ in cohort and preprocessing; coefficient changes cannot be ascribed solely to mutual adjustment. Selected-gene correlations and figure trends are descriptive. No new statistical calculations were performed for this documentation overlay.

## Bounded validation

V2 assessed fixed-input single-patient influence; its supplement added the three non-I primary Stage contrasts. V3 assessed historical QC addbacks with recalculated preprocessing. Amended V4 integrates these and prior complete-case sensitivity evidence, retaining primary versus conditional estimands. The frozen mandatory scope is complete and scientific interpretation QUALIFIED. See [validation](VALIDATION_REPORT.md) and [release scope](REPRODUCIBILITY_AND_RELEASE.md).

Methodological references are listed in [REFERENCES.md](REFERENCES.md). Exact executed parameters and object provenance belong to the frozen v09 code and archived execution records; this narrative is not a replacement for them.
