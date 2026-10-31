# Decision Log

This indexed English record covers all 62 distinct numbered decision labels identified in the historical notebook. It is a source-grounded decision history with current dispositions, not a reconstructed chronology from timestamps. Full surviving English paragraphs and variants are preserved in the separate editorial compendium. Missing labels are D19, D23, D47, D53 and D65; none is invented.

## D01

PDF pages: 96-96.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D01 — Primary Tumor restriction Only RNA-seq files annotated as Primary Tumor were considered.

## D02

PDF pages: 96-96.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D02 — One molecular observation per participant The analytical cohort was designed to contain at most one RNA-seq observation per TCGA participant to avoid pseudoreplication.

## D03

PDF pages: 96-96.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D03 — Multiple Primary Tumor specimens When multiple Primary Tumor specimens existed for one participant, the non-FFPE specimen was retained over the FFPE specimen. Rationale: The distinction was supported by biospecimen FFPE annotations and substantially higher RNA integrity metrics in the non-FFPE specimens.

## D04

PDF pages: 96-97.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D04 — Multiple aliquot-level files Multiple STAR-count files derived from the same RNA analyte were not removed arbitrarily. Targeted expression-level QC was performed before selection.

## D05

PDF pages: 97-97.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D05 — Aliquot selection For the five specimens with two aliquot-level STAR-count files, the aliquot with the greater gene-assignment fraction was retained. The decision was based on expression-level QC rather than aliquot barcode, plate identifier, file size, or alphabetical ordering.

## D06

PDF pages: 97-97.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D06 — CPM sensitivity analysis CPM and log2(CPM + 1) were used only for technical comparison of paired aliquots. This analysis does not define the normalization strategy for downstream statistical modelling.

## D07

PDF pages: 97-98.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D07 — Clinical matching Only participants present in the strictly curated clinical cohort were retained in the matched molecular-clinical cohort. Final matched cohort: N = 1,093.

## D08

PDF pages: 98-98.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D08 — Missing clinical variables Missing age or stage information will be handled separately for each analysis rather than using global complete-case deletion.

## D09

PDF pages: 102-103.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D09 — Canonical RNA-seq manifest The final analytical RNA-seq file set was defined directly from the curated matched molecular-clinical cohort rather than by repeating the original GDC selection criteria independently. This ensures that downstream file acquisition uses exactly the files retained during specimen- and aliquot-level curation. The final manifest contains 1,093 unique participants and 1,093 unique STAR-count files, with one file per participant.

## D10

PDF pages: 104-104.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D10 — Download chunking The final RNA-seq cohort comprised 1,093 files with a total expected size of approximately 4.32 GB. Files were downloaded in moderate-sized chunks rather than as a single large request to improve robustness to network interruptions and simplify recovery from incomplete downloads.

## D11

PDF pages: 106-106.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D11 — Download completeness validation Expression processing was not initiated immediately after download. The downloaded STAR-count file set was first compared against the canonical 1,093-file RNA-seq manifest. Downstream parsing was permitted only after confirming exact file-set completeness, with no missing, unexpected or duplicated files.

## D12

PDF pages: 108-108.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D12 — File-integrity validation All 1,093 downloaded STAR-count files were required to pass MD5 checksum validation against GDC-provided checksums before being used in downstream expression processing. All 1,093 files passed this validation.

## D13

PDF pages: 111-112.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D13 — Gene-level raw count field The `unstranded` STAR gene-count column was retained as the raw count field for construction of the canonical expression matrix. The four STAR summary rows were excluded from the expression matrix. TPM, FPKM and FPKM-UQ values were not used as substitutes for raw counts.

## D14

PDF pages: 112-112.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D14 — Direct matrix construction All 1,093 retained STAR-count files contained an identical ordered set of 60,660 Ensembl gene identifiers. Therefore, raw count vectors were combined directly without gene-ID merging or reordering.

## D15

PDF pages: 119-120.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D15 — Preserve an immutable raw count matrix The 60,660 × 1,093 raw count matrix was retained before expression filtering or normalization. All downstream transformations will be derived from this canonical raw object rather than modifying it in place.

## D16

PDF pages: 120-120.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D16 — Patient-level column identifiers RNA-seq matrix columns were converted from GDC file names to TCGA participant identifiers only after deterministic mapping through the curated RNA-seq manifest. The matrix was then explicitly reordered to match the final clinical metadata table.

## D17

PDF pages: 120-121.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D17 — Globally zero genes Genes with zero counts across all participants were identified but were not removed during raw matrix construction. Gene filtering will be performed as a separate documented preprocessing step using a predefined RNA-seq expression criterion.

## D18

PDF pages: 129-129.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D18 — No automatic sample exclusion from raw-count metrics Initial library size, detected-gene count and zero-fraction metrics were used for sample-level screening only. No sample was excluded solely because it occupied the tail of one raw-count QC distribution. Potential sample exclusions require convergent evidence from subsequent QC analyses.

## D20

PDF pages: 134-135.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D20 — Value-based validation of STAR assigned counts STAR-derived gene-level assigned counts were validated against column sums of the canonical raw count matrix. A difference in R storage mode caused `identical()` to return `FALSE`, but value-based comparison confirmed zero numerical discrepancies across all 1,093 participants.

## D21

PDF pages: 138-139.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D21 — STAR IQR thresholds used for screening, not exclusion Distribution-derived STAR thresholds were used only to prioritize samples for further investigation. STAR categories are not statistically independent, and a single STAR flag was not considered sufficient evidence for sample exclusion.

## D22

PDF pages: 139-139.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D22 — Convergent QC evidence prioritized Samples showing both raw-expression abnormalities and STAR-level technical abnormalities were designated as convergent QC candidates. This designation prioritizes samples for additional investigation and does not itself constitute an exclusion rule.

## D24

PDF pages: 147-147, 192-192.

Historical preliminary filtering snapshot. Final Age and Stage dimensions are 1,072 x 25,551 and 1,063 x 31,157 after the six QC exclusions. Preserve the design-aware filtering rationale, not the preliminary dimensions as final counts.

Historical wording: D24 — Low-expression filtering was performed separately for the Age and Stage analyses because filterByExpr() is design-aware and the two analyses have different eligibility sets and statistical designs.

## D25

PDF pages: 192-192.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D25 — A manually considered CPM ≥1 in ≥20 samples threshold was not adopted. It was retained only as a filtering sensitivity analysis. Design-aware edgeR::filterByExpr() was selected for the primary workflow.

## D26

PDF pages: 192-192.

Historical preliminary filtering snapshot. Final Age and Stage dimensions are 1,072 x 25,551 and 1,063 x 31,157 after the six QC exclusions. Preserve the design-aware filtering rationale, not the preliminary dimensions as final counts.

Historical wording: D26 — Age filtering retained 25,627 genes across 1,078 patients. Stage filtering retained 31,445 genes across 1,069 patients.

## D27

PDF pages: 192-192.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D27 — TMM was selected as the primary between-sample normalization method to account for library-size and compositional differences.

## D28

PDF pages: 193-193.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D28 — TMM normalization factors were treated as compositional adjustments rather than sample-quality scores; extreme factors alone were not used for exclusion.

## D29

PDF pages: 193-193.

limma-voom is the executed primary framework. The proposed edgeR quasi-likelihood comparison was not performed in the frozen V0-V4 protocol.

Historical wording: D29 — limma–voom was selected as the primary inferential framework. edgeR negative-binomial/quasi-likelihood modelling is reserved for later sensitivity/robustness analysis.

## D30

PDF pages: 193-193.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D30 — No separate log2(count + 1) transformation was applied. voom-generated normalized log2-CPM expression values are used for exploratory expression analyses and subsequent limma modelling.

## D31

PDF pages: 193-193.

Historical pre-exclusion checkpoint. Final exclusions and re-estimation occurred later. PCA is diagnostic; its expression-derived evidence is not independent biological replication.

Historical wording: D31 — Post-normalization PCA was performed as an exploratory sample-level QC procedure, not as a formal test of Age or Stage association.

## D32

PDF pages: 193-193.

Historical pre-exclusion checkpoint. Final exclusions and re-estimation occurred later. PCA is diagnostic; its expression-derived evidence is not independent biological replication.

Historical wording: D32 — No sample was excluded solely because of its position in PC1–PC2 space. Final exclusion remains dependent on convergence across raw-expression, STAR, biospecimen, RNA-quality, normalization, and multivariate expression evidence.

## D33

PDF pages: 193-193.

Historical pre-exclusion checkpoint. Final exclusions and re-estimation occurred later. PCA is diagnostic; its expression-derived evidence is not independent biological replication.

Historical wording: D33 — If final sample exclusions occur, filtering, TMM normalization, voom transformation, and final PCA/QC will be rerun before inferential modelling.

## D34

PDF pages: 195-196.

Historical exclusion rationale. The final amended V4 still records incomplete individual exclusion provenance. This narrative and the reported PCA ranks do not remove that limitation; RIN is not the reason for all six exclusions.

Historical wording: D34 — Final exclusion of six technically concordant expression outliers. Six samples previously identified as the most severe STAR-level QC candidates ([participant or aliquot identifier omitted], [participant or aliquot identifier omitted], [participant or aliquot identifier omitted], [participant or aliquot identifier omitted], [participant or aliquot identifier omitted], and [participant or aliquot identifier omitted]) also ranked among the ten most extreme multivariate PCA-distance samples in both the Age and Stage analysis datasets. The concordance was highly consistent across the two analysis-specific expression matrices; [participant or aliquot identifier omitted], for example, ranked first in both analyses. Because these samples showed convergent abnormalities across independent technical and normalized-expression QC layers, all six were excluded from the primary statistical analyses. [participant or aliquot identifier omitted] additionally had independent biospecimen evidence of compromised RNA quality (FFPE, RIN 2.6, 28S/18S ratio 0). Samples displaying unusual PCA behaviour without independent technical abnormalities were retained to avoid removing potentially genuine biological heterogeneity.

## D35

PDF pages: 200-200.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D35 — Final analytical datasets frozen after post-exclusion preprocessing. Following removal of the six QC-supported samples, analysis-specific expression datasets were reconstructed from the canonical raw count matrix and low-expression filtering, TMM normalization, and voom transformation were re-estimated. The final Age dataset contained 1,072 patients and 25,551 retained genes, while the final Stage dataset contained 1,063 patients and 31,157 retained genes. Final PCA showed preservation of broad transcriptomic structure without evidence warranting further technically motivated exclusions. No additional samples will be removed on the basis of exploratory PCA position alone. These datasets were therefore frozen for primary inferential analysis.

## D36

PDF pages: 212-212.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D36—Primary age association model and significance framework. Age at diagnosis was retained as a continuous predictor and tested using limma–voom gene-wise linear modelling. The primary effect estimate was the Age regression coefficient, reported additionally on a 10-year scale for interpretability. Statistical significance was defined using BH-FDR < 0.05. Effect-size thresholds were used descriptively and were not treated as additional hypothesis-testing cutoffs.

## D37

PDF pages: 217-217.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D37 — Transcriptome-wide Age-association visualization. The primary transcriptome-wide visualization used the estimated Age coefficient scaled to a 10-year increase on the x-axis and −log10(BH-FDR) on the y-axis. BH-FDR rather than the unadjusted P-value was displayed to align the visualization with the prespecified multiple-testing criterion used for primary inference. Because Age was modelled as a continuous predictor, the figure was described as an Age-association volcano-style plot rather than a conventional differential-expression volcano plot. Selected labels were used for interpretability and were not intended to constitute a ranking of biological importance.

## D38

PDF pages: 219-219.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D38 — Patient-level visualization of selected Age-associated genes. Four protein-coding genes (ESR1, SOX10, GABRP, and KRT5) were selected to illustrate representative positive and negative Age-associated expression patterns. Patient-level voom log2-CPM expression was plotted against age at diagnosis for all 1,072 patients. Primary effect estimates and BH-FDR values were obtained from the weighted limma-voom model, while Pearson correlations and ordinary least-squares trend lines with 95% confidence intervals were included only as descriptive summaries. Free y-axis scales were used to emphasize within-gene Age-expression patterns; visual slope angles across panels were therefore not interpreted as direct comparisons of effect magnitude.

## D39

PDF pages: 220-220.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D39 — Ranking of Age-associated effect sizes. Among protein-coding genes meeting the primary BH-FDR < 0.05 criterion, genes were ranked by the absolute magnitude of the limma-voom Age coefficient scaled to a 10-year increase in age at diagnosis. The 30 genes with the largest absolute coefficients were visualized. No additional expression-abundance threshold was imposed at this stage because the purpose of the figure was to display the largest statistically supported estimated effects rather than to define a final biological-priority set. Average expression abundance was retained as a separate consideration during downstream interpretation. Twenty-six of the 30 largest absolute effects were negative and four were positive; this imbalance was interpreted specifically as a property of the extreme-effect subset and not as evidence that transcriptome-wide Age associations were predominantly negative.

## D40

PDF pages: 229-229.

This records an earlier interpretation boundary. Later subtype, RIN, purity and combined complete-case sensitivity evidence exists; detailed mechanistic or external validation was not established.

Historical wording: D40 — Scope of preliminary biological interpretation. Biological interpretation within the primary Age analysis was intentionally limited to descriptive contextualization and hypothesis generation. Extensive literature validation, subtype-specific modelling, pathway-level investigation, and mechanistic interpretation were deferred because the primary objective of Phase 1 was to establish a complete and reproducible RNA-seq–clinical association pipeline. These analyses remain candidates for subsequent research-oriented expansion of the project.

## D41

PDF pages: 232-233.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D41 — Final Stage analysis dataset and design. Following removal of the six samples excluded during final expression QC, patients with available pathological stage information were retained for the primary Stage analysis. The resulting dataset contained 1,063 patients: 179 Stage I, 615 Stage II, 249 Stage III, and 20 Stage IV tumors. Stage was represented as a four-level categorical factor with Stage I as the reference level rather than as a numeric predictor, avoiding the assumption of a constant linear expression change across successive stages. Analysis-specific filterByExpr filtering retained 31,157 genes. TMM normalization and voom transformation were re-estimated using the final Stage cohort. The resulting design matrix was full rank, and expression, metadata, and design-matrix sample ordering were identical. Stage IV was retained as a distinct clinical category, although its smaller sample size (n = 20) will be considered when interpreting Stage-IV-specific estimates and contrasts.

## D42

PDF pages: 236-236.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D42 — Primary Stage hypothesis test. The primary Stage association test was defined as a moderated omnibus F-test jointly evaluating all non-reference Stage coefficients. This separates the question of whether expression differs across pathological stages from the subsequent question of which Stage comparisons account for that difference. Signed biological interpretation will therefore be based on explicit contrasts rather than the omnibus statistic.

## D43

PDF pages: 238-238.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D43 — Stage multiple-testing framework. The primary Stage-associated gene set was defined by BH-FDR < 0.05 on the moderated omnibus F-test. No additional effect-size threshold was imposed for statistical significance. More stringent FDR thresholds were summarized descriptively but were not used to redefine the primary gene set. Directional interpretation was deferred to contrast-specific analyses.

## D44

PDF pages: 239-239.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D44 — Pairwise characterization of pathological-stage associations. Following the primary omnibus Stage test, all six pairwise contrasts among Stages I–IV were evaluated to characterize the location, magnitude, and direction of Stage-associated expression differences. BH correction was applied separately across the 31,157 tested genes within each contrast. The omnibus test remained the primary definition of overall Stage association, whereas pairwise contrasts were used for contrast-specific inference and interpretation. The number of significant genes differed substantially across comparisons, with 1,004 genes identified for Stage II versus I, 36 for Stage III versus I, 57 for Stage IV versus I, 646 for Stage III versus II, 5 for Stage IV versus II, and 7 for Stage IV versus III at BH-FDR < 0.05. These differences were not interpreted as direct measures of biological distance between stages because statistical power and precision differed across comparisons, particularly for Stage IV (n = 20).

## D45

PDF pages: 240-241.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D45 — Annotation of Stage-association results. Gene identifiers from the primary Stage analysis were annotated using the GENCODE v36 gene annotation embedded in the original GDC STAR-count files, thereby maintaining consistency with the expression quantification source. All 31,157 tested genes were successfully annotated without missing gene names or biotypes. Among the 1,285 genes meeting the primary omnibus BH-FDR < 0.05 criterion, 972 were protein-coding and 232 were lncRNAs. Biotype-specific proportions were summarized descriptively and were not interpreted as formal enrichment statistics.

## D46

PDF pages: 241-241.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D46 — Separation of statistical significance, effect magnitude, and expression abundance. Omnibus BH-FDR defined statistical evidence for overall Stage association, whereas the maximum absolute Stage-versus-I coefficient was used to describe effect magnitude. Average expression abundance was considered separately because large coefficients could occur for low-abundance genes. No post-hoc abundance or effect-size threshold was added to the primary significance definition.

## D48

PDF pages: 243-243.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D48 — Patient-level visualization of selected Stage-associated genes. Four protein-coding genes (ACKR1, ADH1B, ADIPOQ, and PIR) were selected to illustrate representative Stage-associated expression patterns among the 1,063 tumors. Patient-level voom log2-CPM expression distributions were displayed across the four pathological-stage groups using boxplots overlaid with individual tumor observations. Omnibus BH-FDR values were reported from the primary moderated Stage test. The plots were descriptive visualizations of the fitted Stage associations rather than independent statistical tests. Free y-axis scales were used to improve within-gene visualization and therefore visual differences in vertical scale across gene panels were not interpreted as comparisons of effect magnitude. The substantially smaller Stage IV group (n = 20) was retained and explicitly reported.

## D49

PDF pages: 243-244.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D49 — Visualization of Stage-specific coefficient patterns. To characterize the structure of pathological-stage associations, the 30 protein-coding genes with the strongest omnibus statistical evidence were visualized using their fitted Stage II-versus-I, Stage III-versus-I, and Stage IV-versus-I limma-voom coefficients. Raw model coefficients were displayed without row-wise standardization so that both direction and relative effect magnitude were retained. Genes were ordered by omnibus BH-FDR rather than clustered, and Stage comparisons were maintained in their clinical order. For visualization, the diverging color scale was centered at zero and capped symmetrically at ±2.5 voom log2-CPM; coefficients outside this range were color-saturated but were not altered in the underlying data.

## D50

PDF pages: 250-250.

This records an earlier interpretation boundary. Later subtype, RIN, purity and combined complete-case sensitivity evidence exists; detailed mechanistic or external validation was not established.

Historical wording: D50 — Scope of preliminary Stage interpretation. Biological interpretation of the primary Stage analysis was intentionally restricted to descriptive pattern characterization and hypothesis generation. The analysis established transcriptomic heterogeneity associated with pathological stage but did not attempt to infer causal progression mechanisms. Detailed literature validation, molecular-subtype adjustment, tumor-composition modelling, pathway analysis, and ordinal Stage-trend testing were deferred to subsequent research-oriented extensions of the project.

## D51

PDF pages: 264-265.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D51 — Definition of the common Age–Stage gene universe. Integration of the primary Age and Stage analyses was restricted to genes tested in both analyses to ensure comparable statistical opportunity for association with either clinical variable. All 25,551 genes tested in the Age analysis were also present among the 31,157 genes tested in the Stage analysis, whereas 5,606 additional genes were tested only in the Stage analysis because of analysis-specific low-expression filtering. The common 25,551-gene set was therefore defined as the reference universe for direct Age–Stage comparisons, overlap analyses, and classification of genes as Age-associated, Stage-associated, associated with both, or associated with neither variable.

## D52

PDF pages: 265-265.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D52 — Construction of the integrated Age–Stage results table. Primary Age and Stage results were aligned using the 25,551 genes tested in both analyses. Gene identifiers, gene symbols, and gene biotypes were fully concordant between the two result sets, with no duplicated identifiers or missing Age or Stage FDR values. Within this common universe, 9,063 genes met the primary Age BH-FDR < 0.05 criterion and 1,166 met the primary Stage omnibus BH-FDR < 0.05 criterion. Of the 1,285 Stage-associated genes identified in the full Stage analysis, 119 were outside the common universe because they had not passed the analysis-specific filtering criterion for the Age analysis and were therefore excluded from direct Age–Stage comparison.

## D54

PDF pages: 267-267.

Historical descriptive overlap inference. The Fisher estimate and P value are retained as source wording, not recomputed here. Correlated genes and post-selection limit interpretation; this is not a new confirmatory conjunction-FDR result.

Historical wording: D54 — Enrichment of overlap between Age- and Stage-associated gene sets. Among the 25,551 genes tested in both primary analyses, 513 were significant for both Age and Stage, compared with approximately 413.6 genes expected under independence, corresponding to an observed-to-expected overlap ratio of 1.24. Membership in the Age-associated and Stage-associated gene sets was positively associated (Fisher's exact-test odds ratio = 1.46, 95% CI 1.29–1.64; P = 9.07 × 10 − 1 0 ). This result was interpreted as evidence of greater-than-expected statistical overlap between the two association sets, rather than as evidence of a shared causal mechanism. Because transcriptomic features are correlated and the analysis dichotomizes continuous statistical evidence at BH-FDR < 0.05, the overlap test was treated as a descriptive inferential summary rather than a mechanistic test.

## D55

PDF pages: 269-269.

Historical descriptive characterization. Verify individual coefficients, biotypes and summary values against saved tables before using them in a new numerical claim.

Historical wording: D55 — Characterization of genes associated with both Age and Stage. Among the 513 genes meeting BH-FDR < 0.05 in both primary analyses, 339 (66.1%) had negative Age coefficients and 174 (33.9%) had positive Age coefficients. Stage-specific coefficients showed heterogeneous patterns: 141 genes (27.5%) had negative coefficients for all three Stage-versus-I comparisons, 191 (37.2%) had positive coefficients for all three comparisons, and 181 (35.3%) showed mixed coefficient signs across Stage comparisons. Thus, statistical overlap between the Age and Stage association sets did not correspond to a single concordant direction of expression change. The shared set consisted predominantly of protein-coding genes (440/513; 85.8%), with lncRNAs accounting for 60/513 (11.7%); these proportions were treated descriptively rather than as evidence of biotype enrichment.

## D56

PDF pages: 269-269.

Historical descriptive characterization. Verify individual coefficients, biotypes and summary values against saved tables before using them in a new numerical claim.

Historical wording: D56 — Effect magnitude and biotype representation among shared Age–Stage associations. Among the 513 genes associated with both Age and Stage, the median absolute Age coefficient was 0.103 voom log2-CPM per 10 years (IQR 0.073–0.152), while the median maximum absolute Stage-versus-I coefficient was 0.489 log2-CPM (IQR 0.347–0.665). Protein-coding genes accounted for 440 shared associations and represented 2.63% of all protein-coding genes tested in the common universe, whereas 60 lncRNAs were shared, representing 1.05% of tested lncRNAs. These denominator-aware proportions were treated descriptively and were not interpreted as formal biotype enrichment.

## D57

PDF pages: 270-271.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D57 — Definition of the joint Age–Stage analysis cohort. The secondary joint analysis was restricted to patients represented in both final primary analysis cohorts after the six QC-supported exclusions. Patient matching was performed using the canonical TCGA patient_id rather than data-frame row names. The final Age cohort contained 1,072 unique patients and the final Stage cohort contained 1,063 unique patients, with 1,048 patients shared between them. Twenty-four patients were present only in the Age cohort and 15 only in the Stage cohort. All 1,048 joint-cohort patient identifiers were uniquely represented in the final raw count matrix, and Age and Stage clinical variables were concordant across the two source metadata tables. The joint cohort was therefore frozen at 1,048 patients for secondary mutually adjusted modelling.

## D58

PDF pages: 271-272.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D58 — Specification of the joint Age–Stage model. The joint analysis included 1,048 patients with complete Age and Stage information following the final QC exclusions. The cohort comprised 178 Stage I, 603 Stage II, 247 Stage III, and 20 Stage IV tumors. Age at diagnosis ranged from 26.57 to 90.00 years and was parameterized in 10-year units for direct coefficient interpretation. Pathological stage was modelled as a four-level categorical factor with Stage I as the reference level. The joint design matrix contained five coefficients—intercept, Age per 10 years, and Stage II, III, and IV versus Stage I—and was full rank (rank 5/5). Age distributions overlapped across Stage groups, with mean ages ranging from 58.69 to 60.31 years. The joint model was defined as a secondary mutually adjusted analysis and did not replace the frozen primary Age-only and Stage-only analyses.

## D59

PDF pages: 272-272.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D59 — Joint Age–Stage expression preprocessing. The joint analysis was reconstructed from the final canonical raw count matrix rather than reusing expression values or precision weights from either primary analysis. Raw counts for the 1,048-patient joint cohort comprised 60,660 genes, with exact patient-order concordance between the count matrix, metadata, and joint design matrix. Joint-design-specific filterByExpr filtering retained 31,285 genes and removed 29,375. TMM normalization was then re-estimated within the joint cohort, followed by joint-design-specific voom transformation. The resulting expression and precision-weight matrices contained 31,285 genes across 1,048 patients, with no missing values and exact alignment to the design matrix. The voom mean–variance trend was considered suitable for downstream limma modelling. The joint-specific preprocessing outputs were therefore frozen for secondary mutually adjusted inference.

## D60

PDF pages: 273-273.

Retained historical documentation wording. Reuse only with the current v09 model, cohort, unit and testing-family definitions. The current project-wide validation decision is COMPLETE within the frozen mandatory scope, with scientific interpretation QUALIFIED.

Historical wording: D60 — Mutually adjusted Age and Stage inference. The secondary joint limma–voom model was fitted to 31,285 genes across 1,048 patients using Age per 10 years and pathological Stage as simultaneous predictors. The model had five coefficients and 1,043 residual degrees of freedom for every tested gene. After adjustment for Stage, 9,413 genes met BH-FDR < 0.05 for Age, including 5,470 negative and 3,943 positive Age coefficients. After adjustment for Age, 1,201 genes met BH-FDR < 0.05 in the three-degree-of-freedom omnibus Stage test. These discovery counts were not directly compared with those from the primary Age-only and Stage-only analyses because the joint analysis used a different patient cohort and gene-testing universe. Primary-versus-adjusted stability was therefore evaluated separately on matched gene sets.

## D61

PDF pages: 274-274.

Concordance under primary versus Joint specifications is a limited descriptive comparison. Cohort and preprocessing also differ. It does not isolate adjustment, establish equivalence, or replace final influence, addback and subtype qualifications.

Historical wording: D61 — Stability of Age associations after adjustment for pathological Stage. Primary and Stage-adjusted Age associations were compared across the same 25,551 genes tested in the primary Age analysis. Age effect estimates were highly concordant between the primary and adjusted models (Pearson r = 0.9979; Spearman ρ = 0.9968), with coefficient direction preserved for 96.33% of genes. The median absolute change in the Age coefficient was 0.0049 voom log2-CPM per 10 years, and 99% of genes changed by less than 0.0223; 99.99% changed by less than 0.05. Of the 9,063 genes significant in the primary Age analysis, 8,038 (88.69%) remained significant at BH-FDR < 0.05 after Stage adjustment; 1,025 were significant only in the primary analysis, while 301 genes that were not significant in the primary analysis became significant after adjustment. Together, these results indicate that Age effect estimates were highly stable to adjustment for pathological Stage in this cohort. Changes in significance status were interpreted cautiously because the adjusted analysis also differed in cohort size and preprocessing and because statistical significance depends on both effect magnitude and precision.

## D62

PDF pages: 275-276.

Concordance under primary versus Joint specifications is a limited descriptive comparison. Cohort and preprocessing also differ. It does not isolate adjustment, establish equivalence, or replace final influence, addback and subtype qualifications.

Historical wording: D62 — Stability of Stage associations after adjustment for Age. Primary and Age-adjusted Stage associations were compared across 31,154 genes tested in both analyses. Of the 1,285 genes significant in the primary omnibus Stage analysis, 1,117 (86.93%) remained significant at BH-FDR < 0.05 after adjustment for Age; 168 were significant only in the primary analysis, while 84 genes became significant only after adjustment. Stage-specific coefficient estimates were highly concordant between the primary and adjusted models: Pearson correlations were 0.9949, 0.9961, and 0.9994 for Stage II–I, III–I, and IV–I, respectively, with corresponding sign concordance of 96.32%, 96.97%, and 98.60%. Median absolute coefficient changes were 0.0106, 0.0092, and 0.0106 voom log2-CPM, respectively. These results indicate that Stage effect estimates were highly stable to adjustment for Age. As with the Age stability analysis, changes in statistical significance were interpreted cautiously because the adjusted analysis differed modestly in cohort composition and preprocessing, and threshold crossing does not by itself imply that Age explains or removes a Stage association.

## D63

PDF pages: 278-278.

Concordance under primary versus Joint specifications is a limited descriptive comparison. Cohort and preprocessing also differ. It does not isolate adjustment, establish equivalence, or replace final influence, addback and subtype qualifications.

Historical wording: D63 — Stability of shared Age–Stage associations under mutual adjustment. The 513 genes identified as significant for both Age and Stage in the primary analyses were followed into the secondary joint model. All 513 genes were represented in the joint testing universe. After mutual adjustment, 418 genes (81.48%) remained significant for both Age and Stage at BH-FDR < 0.05. An additional 63 genes remained significant for Age but not Stage, 25 remained significant for Stage but not Age, and only 7 genes were no longer significant for either variable. Overall, 481/513 (93.76%) retained an Age association and 443/513 (86.35%) retained a Stage association in the joint model. These results indicate substantial statistical stability of the primary shared association set under mutual adjustment. Because the 513 genes were selected on the basis of significance in the primary analyses, this follow-up was treated as descriptive characterization rather than an independent significance test.

## D64

PDF pages: 280-280.

Use the archived 25,551-gene integrated master and Figure 8 metrics. Do not mix them with the 31,154-gene Stage comparison mentioned elsewhere. Strong global correlation cannot certify every candidate contrast.

Historical wording: D64 — Construction and validation of the final integrated results table. A canonical integrated results table was constructed for the 25,551 genes tested in both primary Age and Stage analyses. The table contained 34 columns covering gene annotation, primary Age and Stage statistics, mutually adjusted Age and Stage statistics, primary and adjusted association categories, and category transitions. Gene identifiers were unique, non-missing, and exactly aligned across the primary and joint results. Primary classification reproduced 8,550 Age-only, 653 Stage-only, 513 shared, and 15,835 neither-associated genes. Within the common universe, the adjusted analysis identified 8,339 Age-associated and 1,083 Stage-associated genes, including 479 associated with both variables. Of the 513 genes associated with both variables in the primary analyses, 418 retained both associations after mutual adjustment, reproducing the previously validated shared-gene stability results. The integrated table was designated as the canonical source for downstream integrated visualization and reporting.

## D66

PDF pages: 281-281.

Use the archived 25,551-gene integrated master and Figure 8 metrics. Do not mix them with the 31,154-gene Stage comparison mentioned elsewhere. Strong global correlation cannot certify every candidate contrast.

Historical wording: D66 — Primary versus Adjusted Effect-Estimate Stability Visualization Decision Log · Validated Figure 8 visualized gene-level coefficient stability between primary and mutually adjusted Age–Stage models across the common universe of 25,551 genes. Four panels represented Age per decade and Stage II, III, and IV relative to Stage I. Pearson correlations ranged from 0.99485 to 0.99936, with sign concordance between 96.29% and 98.57%. Median absolute coefficient differences ranged from 0.00490 to 0.01092. The visualizations showed strong concordance of effect estimates across model specifications. However, comparisons also involved changes in cohort composition and preprocessing, and therefore do not isolate the effect of covariate adjustment or establish causal independence.

## D67

PDF pages: 281-282.

Historical Figure 9 selection definition. The 30-gene panel is descriptive prioritization from 418 shared retained associations. Final scientific interpretation remains QUALIFIED; within-column colors are not comparable between contrasts.

Historical wording: D67 — Integrated Age–Stage Coefficient Heatmap Decision Log · Validated Figure 9 was constructed using 30 genes selected from 418 stable shared Age–Stage associations. Genes were prioritized through a deterministic worst-rank procedure combining adjusted Age and Stage FDR rankings, with rank-sum and Ensembl ID tie-breaking. The heatmap displayed mutually adjusted Age coefficients per decade and Stage II, III, and IV coefficients relative to Stage I. For visualization, coefficients were scaled independently within each contrast by the maximum absolute coefficient among the 30 selected genes, preserving coefficient signs and zero. Gene ordering followed the predefined selection rank without clustering. The figure was visually reviewed and validated. Color intensities represent within-contrast scaled coefficients and are not directly comparable across contrasts.
