# Study design and clinical curation

## Research questions and analytical unit

This project examines cross-sectional associations between bulk tumor gene expression, age at diagnosis and pathological stage in TCGA-BRCA. The analytical unit is the participant. One curated Primary Tumor RNA-seq file is retained per participant to avoid treating related specimens or aliquots as independent observations.

The primary Age model tests a continuous age coefficient. The primary Stage model tests expression heterogeneity across Stages I-IV. A secondary Joint model includes both age and stage. The Joint analysis complements the separate primary analyses; it does not replace them.

| Analysis | Question | Interpretation boundary |
|---|---|---|
| Age | Is expression associated with age at diagnosis? | Between-tumor association, not within-person ageing |
| Stage | Does expression differ across pathological-stage groups? | Overall heterogeneity and explicitly defined contrasts |
| Joint | Does an association remain under mutual Age-Stage adjustment? | Adjustment for these two variables, not all confounders |

Age was derived from diagnosis age in days divided by 365.25 and retained continuously. The primary fit uses years; effects are also reported per decade. Joint-model age is parameterized in decades. The unit must accompany every coefficient.

Stage is a four-level unordered factor with Stage I as the reference. This does not assume equal numerical differences between successive stages. Earlier notebook proposals for numeric or ordinal encoding are historical planning material.

Source notebook: PDF pages 27-45, 56-59 and 208. Current evidence: E01-E04. Proposed mechanistic, longitudinal or clinical-prediction interpretations are outside the demonstrated scope.



## Clinical curation and stage harmonization

The archived clinical table contained 5,529 rows and 202 fields representing 1,095 unique participants. Repeated clinical rows were not treated as independent patients. The retained cohort required an explicitly designated primary-disease diagnosis. This produced 1,093 participants with one unique primary diagnosis each.

Substages were mapped to their corresponding major stages: I/IA/IB to I; II/IIA/IIB to II; III/IIIA/IIIB/IIIC to III; and IV to IV. The original value was retained separately from the derived stage factor. Stage X was considered unavailable for stage-specific analysis rather than a fifth biological stage.

| Clinical stage of processing | Participants | Meaning |
|---|---:|---|
| Unique participants in clinical export | 1,095 | Before strict diagnosis selection |
| Curated clinical and matched molecular cohort | 1,093 | Starting cohort before expression QC exclusions |
| Age available before QC exclusions | 1,078 | 15 participants lacked usable age |
| Interpretable stage before QC exclusions | 1,069 | 13 Stage X and 11 missing-stage records |
| Both age and interpretable stage before QC | 1,054 | Eligibility before the six QC exclusions |

The pre-QC stage distribution was 181/619/249/20 for Stages I/II/III/IV. These counts describe an earlier processing stage. The final primary Stage distribution is 179/615/249/20; the Joint distribution is 178/603/247/20.

Clinical records and individual participant identifiers are not distributed in this documentation package. The historical vignette concerning an explicit-primary-diagnosis exception is summarized without identifying details in the editorial compendium.

Source notebook: PDF pages 48-62 and 89-95. Current evidence: E01-E03. Historical clinical curation is distinct from subsequent QC and analysis-specific cohort selection.



## Analysis-specific eligibility and final cohorts

Missingness was handled according to the variables required by each analysis. A participant lacking stage could still contribute to the Age analysis. Globally restricting every analysis to complete cases would change its population unnecessarily.

| Model | Before QC exclusions | Final patients | Final tested genes |
|---|---:|---:|---:|
| Age | 1,078 | 1,072 | 25,551 |
| Stage | 1,069 | 1,063 | 31,157 |
| Joint Age-Stage | 1,054 | 1,048 | 31,285 |

Six historical QC exclusions were applied before definitive primary modelling. Their shared rationale was concordant technical and multivariate expression-space abnormalities. The final validation retains incomplete individual exclusion provenance; the notebook narrative is not a substitute for missing contemporaneous evidence.

The primary Age model was not a sex- or batch-adjusted analysis, and the primary Stage model did not simultaneously adjust for Age. Later RIN, purity, subtype and combined complete-case models answer additional questions with their own eligible cohorts. They do not retroactively change the primary model definitions.

The primary and Joint models differ in cohort composition and preprocessing as well as formula. Changes between their coefficients or significance labels therefore cannot be attributed solely to covariate adjustment. Same-patient comparisons are reported separately where available.

An earlier notebook branch suggested retaining a participant without an explicit primary-disease flag; the consolidated curation report and frozen implementation follow the strict rule. This discrepancy is retained in the editorial history rather than combined into the current Methods.

Source notebook: PDF pages 49, 55, 197-200 and 270-273. Current evidence: E01-E03 and E10. These are observational bulk-expression associations; no causal mechanism, within-patient progression or clinical utility is established.
