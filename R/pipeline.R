# Phase 3.3 v09: real-data-tested callable analysis; sourcing alone performs no analysis.
# Use run.R with explicit --execute and separate input/output roots.
# Full real-data reproducibility remains Phase 3.3.
tcga_brca_analysis <- function(ctx) {
analysis_env <- environment()
read.delim <- function(file, ...) pipeline_read(ctx, file, utils::read.delim, ...)
read.csv <- function(file, ...) pipeline_read(ctx, file, utils::read.csv, ...)
readRDS <- function(file, ...) pipeline_read(ctx, file, base::readRDS, ...)
write.csv <- function(x, file, ...) pipeline_publish(ctx, file,
  function(p) utils::write.csv(x, file = p, ...),
  producer = paste(deparse(substitute(x)), collapse = " "),
  settings = list(format = "CSV", dimensions = dim(x), arguments = list(...)))
saveRDS <- function(object, file, ...) pipeline_publish(ctx, file,
  function(p) {
    has_plots <- is.list(object) && any(vapply(object, inherits, logical(1), what = "ggplot"))
    portable <- if (has_plots) pipeline_detach_analysis(object, analysis_env, bindings = "figure7_group_order") else object
    base::saveRDS(portable, file = p, ...)
  },
  producer = paste(deparse(substitute(object)), collapse = " "),
  settings = list(format = "RDS", class = class(object), dimensions = dim(object), arguments = list(...)))
writeLines <- function(text, con, ...) pipeline_publish(ctx, con,
  function(p) base::writeLines(text, con = p, ...),
  producer = paste(deparse(substitute(text)), collapse = " "), settings = list(format = "text"))
ggsave <- function(filename, plot, ...) pipeline_publish(ctx, filename,
  function(p) ggplot2::ggsave(filename = p, plot = plot,
    device = if (tolower(tools::file_ext(filename)) == "pdf") pipeline_pdf_device else NULL, ...),
  producer = paste(deparse(substitute(plot)), collapse = " "),
  settings = list(format = tools::file_ext(filename), arguments = list(...), RNG_seed = ctx$config$seed))

library(TCGAbiolinks)
input_dir <- ctx$input_dir
output_dir <- ctx$output_dir
s9_export_dir <- file.path(output_dir, "data_processed", "S9_integrated_results")
clinical_dir <- file.path(input_dir, "data_raw", "clinical")
pipeline_event(ctx, "clinical")

clinical_raw <- read.delim(
  file.path(clinical_dir, "clinical.tsv"),
  check.names = FALSE,
  na.strings = c("'--", "")
)

dim(clinical_raw)
names(clinical_raw)
grep(
  "case|submitter|patient",
  names(clinical_raw),
  value = TRUE,
  ignore.case = TRUE
)
grep(
  "age",
  names(clinical_raw),
  value = TRUE,
  ignore.case = TRUE
)
grep(
  "sex|gender",
  names(clinical_raw),
  value = TRUE,
  ignore.case = TRUE
)
grep(
  "primary",
  names(clinical_raw),
  value = TRUE,
  ignore.case = TRUE
)
# 1. Number of rows and unique patients
nrow(clinical_raw)

length(unique(clinical_raw$cases.case_id))
length(unique(clinical_raw$cases.submitter_id))


# 2. Confirm project(s)
table(clinical_raw$project.project_id, useNA = "ifany")


# 3. How many clinical rows does each patient have?
rows_per_patient <- table(clinical_raw$cases.submitter_id)

summary(as.numeric(rows_per_patient))

table(rows_per_patient)


# 4. Inspect the primary-diagnosis indicator
table(
  clinical_raw$diagnoses.diagnosis_is_primary_disease,
  useNA = "ifany"
)


# 5. Inspect target clinical variables
summary(clinical_raw$diagnoses.age_at_diagnosis)

table(
  clinical_raw$diagnoses.ajcc_pathologic_stage,
  useNA = "ifany"
)

table(
  clinical_raw$demographic.sex_at_birth,
  useNA = "ifany"
)
# Number of unique diagnosis IDs
length(unique(clinical_raw$diagnoses.diagnosis_id))
# Number of unique primary diagnosis IDs
primary_rows <- clinical_raw[
  clinical_raw$diagnoses.diagnosis_is_primary_disease == "true" &
    !is.na(clinical_raw$diagnoses.diagnosis_is_primary_disease),
]

length(unique(primary_rows$diagnoses.diagnosis_id))
length(unique(primary_rows$cases.submitter_id))
# How many unique primary diagnoses does each patient have?
primary_diag_per_patient <- tapply(
  primary_rows$diagnoses.diagnosis_id,
  primary_rows$cases.submitter_id,
  function(x) length(unique(x))
)

table(primary_diag_per_patient)
# Are age and stage consistent across repeated rows
# belonging to the same primary diagnosis?

age_values_per_diag <- tapply(
  primary_rows$diagnoses.age_at_diagnosis,
  primary_rows$diagnoses.diagnosis_id,
  function(x) length(unique(na.omit(x)))
)

stage_values_per_diag <- tapply(
  primary_rows$diagnoses.ajcc_pathologic_stage,
  primary_rows$diagnoses.diagnosis_id,
  function(x) length(unique(na.omit(x)))
)

table(age_values_per_diag)
table(stage_values_per_diag)
# Identify patients without a TRUE primary-disease diagnosis

all_patients <- unique(clinical_raw$cases.submitter_id)

primary_patients <- unique(primary_rows$cases.submitter_id)

patients_without_primary <- setdiff(
  all_patients,
  primary_patients
)

patients_without_primary
clinical_raw[
  clinical_raw$cases.submitter_id %in% patients_without_primary,
  c(
    "cases.submitter_id",
    "diagnoses.diagnosis_id",
    "diagnoses.diagnosis_is_primary_disease",
    "diagnoses.primary_diagnosis",
    "diagnoses.age_at_diagnosis",
    "diagnoses.ajcc_pathologic_stage"
  )
]
unique(
  clinical_raw[
    clinical_raw$cases.submitter_id %in% patients_without_primary,
    c(
      "cases.submitter_id",
      "diagnoses.diagnosis_id",
      "diagnoses.diagnosis_is_primary_disease",
      "diagnoses.primary_diagnosis",
      "diagnoses.age_at_diagnosis",
      "diagnoses.ajcc_pathologic_stage"
    )
  ]
)
primary_clinical <- unique(
  primary_rows[
    ,
    c(
      "cases.case_id",
      "cases.submitter_id",
      "diagnoses.diagnosis_id",
      "diagnoses.age_at_diagnosis",
      "diagnoses.ajcc_pathologic_stage",
      "demographic.sex_at_birth"
    )
  ]
)

dim(primary_clinical)
# Rename for clarity
clinical_patient <- primary_clinical

# Check that each patient appears exactly once
nrow(clinical_patient)

length(unique(clinical_patient$cases.submitter_id))

anyDuplicated(clinical_patient$cases.submitter_id)
clinical_patient$age_years <-
  clinical_patient$diagnoses.age_at_diagnosis / 365.25

summary(clinical_patient$age_years)
sum(is.na(clinical_patient$diagnoses.age_at_diagnosis))

sum(is.na(
  clinical_patient$diagnoses.ajcc_pathologic_stage
))
sum(
  !is.na(clinical_patient$age_years) &
    !is.na(clinical_patient$diagnoses.ajcc_pathologic_stage)
)
table(
  clinical_patient$diagnoses.ajcc_pathologic_stage,
  useNA = "ifany"
)
# Keep the original pathological stage untouched
clinical_patient$stage_raw <-
  clinical_patient$diagnoses.ajcc_pathologic_stage


# Create harmonized major stage groups
clinical_patient$stage_group <- NA_character_

clinical_patient$stage_group[
  clinical_patient$stage_raw %in%
    c("Stage I", "Stage IA", "Stage IB")
] <- "Stage I"

clinical_patient$stage_group[
  clinical_patient$stage_raw %in%
    c("Stage II", "Stage IIA", "Stage IIB")
] <- "Stage II"

clinical_patient$stage_group[
  clinical_patient$stage_raw %in%
    c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")
] <- "Stage III"

clinical_patient$stage_group[
  clinical_patient$stage_raw == "Stage IV"
] <- "Stage IV"
table(
  clinical_patient$stage_group,
  useNA = "ifany"
)
table(
  clinical_patient$stage_raw,
  clinical_patient$stage_group,
  useNA = "ifany"
)
# ------------------------------------------------------------
# Define analysis-specific eligibility flags
# ------------------------------------------------------------

# Eligible for age-based analyses
clinical_patient$age_available <-
  !is.na(clinical_patient$age_years)

# Eligible for pathological-stage analyses
clinical_patient$stage_available <-
  !is.na(clinical_patient$stage_group)

# Eligible for analyses requiring both variables
clinical_patient$age_stage_available <-
  clinical_patient$age_available &
  clinical_patient$stage_available
table(clinical_patient$age_available)
table(clinical_patient$stage_available)
table(clinical_patient$age_stage_available)
# Cross-tabulation of analysis eligibility
table(
  Age = clinical_patient$age_available,
  Stage = clinical_patient$stage_available
)
# ============================================================
# Create analysis-ready patient-level clinical metadata
# ============================================================

clinical_clean <- data.frame(
  
  # GDC identifiers
  case_id =
    clinical_patient$cases.case_id,
  
  patient_id =
    clinical_patient$cases.submitter_id,
  
  diagnosis_id =
    clinical_patient$diagnoses.diagnosis_id,
  
  # Age at diagnosis
  age_at_diagnosis_days =
    clinical_patient$diagnoses.age_at_diagnosis,
  
  age_years =
    clinical_patient$age_years,
  
  # Pathological stage
  stage_raw =
    clinical_patient$stage_raw,
  
  stage_group =
    clinical_patient$stage_group,
  
  # Demographic information
  sex =
    clinical_patient$demographic.sex_at_birth,
  
  # Analysis eligibility flags
  age_available =
    clinical_patient$age_available,
  
  stage_available =
    clinical_patient$stage_available,
  
  age_stage_available =
    clinical_patient$age_stage_available,
  
  stringsAsFactors = FALSE
)
# ============================================================
# Final integrity checks
# ============================================================

# Dimensions
dim(clinical_clean)

# Patient IDs must be unique
anyDuplicated(clinical_clean$patient_id)

# Patient IDs must not be missing
sum(is.na(clinical_clean$patient_id))

# Age summary
summary(clinical_clean$age_years)

# Harmonized pathological stage
table(
  clinical_clean$stage_group,
  useNA = "ifany"
)

# Sex distribution
table(
  clinical_clean$sex,
  useNA = "ifany"
)

# Eligibility summary
c(
  total_curated = nrow(clinical_clean),
  age_eligible = sum(clinical_clean$age_available),
  stage_eligible = sum(clinical_clean$stage_available),
  age_stage_eligible =
    sum(clinical_clean$age_stage_available)
)

# Inspect final structure
str(clinical_clean)
# Set pathological stage as an unordered categorical factor
clinical_clean$stage_group <- factor(
  clinical_clean$stage_group,
  levels = c("Stage I", "Stage II", "Stage III", "Stage IV"),
  ordered = FALSE
)

# Optional: sex as categorical factor
clinical_clean$sex <- factor(
  clinical_clean$sex
)

# Confirm final structure
str(clinical_clean)

# Confirm stage levels
levels(clinical_clean$stage_group)

# Final stage count check
table(clinical_clean$stage_group, useNA = "ifany")
# ============================================================
# Save curated clinical metadata
# ============================================================

processed_dir <- file.path(
  output_dir,
  "data_processed"
)

dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

pipeline_assert_ids(clinical_clean$patient_id, "curated clinical patients")
stopifnot(all(clinical_raw$project.project_id == "TCGA-BRCA"))
# Human-readable version
write.csv(
  clinical_clean,
  file = file.path(
    processed_dir,
    "TCGA_BRCA_clinical_curated.csv"
  ),
  row.names = FALSE,
  na = ""
)

# R-native version preserving data types
saveRDS(
  clinical_clean,
  file = file.path(
    processed_dir,
    "TCGA_BRCA_clinical_curated.rds"
  )
)
# ============================================================
# Validate saved RDS file
# ============================================================

clinical_check <- readRDS(
  file.path(
    processed_dir,
    "TCGA_BRCA_clinical_curated.rds"
  )
)

stopifnot(identical(clinical_clean, clinical_check))
# ============================================================
pipeline_event(ctx, "cohort_selection")
# STEP 2 — RNA-seq cohort and file metadata curation
# Step 2A: Query TCGA-BRCA STAR-Counts Primary Tumor data
# ============================================================

query_exp <- pipeline_query(ctx, "query_exp", list(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor",
  access = "open"
), kind = "query")
# Obtain the file/sample metadata returned by the GDC query
exp_results <- getResults(query_exp)
stopifnot(all(c("file_id", "file_name", "cases.submitter_id", "sample.submitter_id",
                "cases", "file_size", "sample_type", "is_ffpe") %in% names(exp_results)))
pipeline_assert_ids(exp_results$file_id, "GDC file IDs")

dim(exp_results)

names(exp_results)
# Number of files returned
nrow(exp_results)

# Inspect the first few records
head(exp_results)
grep(
  "case|submitter|sample|barcode|file|portion|aliquot",
  names(exp_results),
  value = TRUE,
  ignore.case = TRUE
)
# ============================================================
# STEP 2B — Audit RNA-seq file, sample, and patient structure
# ============================================================

# Total STAR-Counts files
nrow(exp_results)

# Unique GDC file UUIDs
length(unique(exp_results$file_id))

# Unique detailed TCGA barcodes
length(unique(exp_results$cases))

# Unique sample/vial barcodes
length(unique(exp_results$sample.submitter_id))

# Unique patients
length(unique(exp_results$cases.submitter_id))
# File IDs should be unique
anyDuplicated(exp_results$file_id)

# File names should also be unique
anyDuplicated(exp_results$file_name)

# Are sample/vial barcodes duplicated?
anyDuplicated(exp_results$sample.submitter_id)

# Are patients represented more than once?
anyDuplicated(exp_results$cases.submitter_id)
# Number of RNA-seq files per patient
files_per_patient <- table(
  exp_results$cases.submitter_id
)

summary(as.numeric(files_per_patient))

table(files_per_patient)
# Number of unique sample/vial barcodes per patient
samples_per_patient <- tapply(
  exp_results$sample.submitter_id,
  exp_results$cases.submitter_id,
  function(x) length(unique(x))
)

summary(samples_per_patient)

table(samples_per_patient)
patients_multiple_samples <- names(
  samples_per_patient[
    samples_per_patient > 1
  ]
)

length(patients_multiple_samples)

patients_multiple_samples
duplicate_exp_metadata <- exp_results[
  exp_results$cases.submitter_id %in%
    patients_multiple_samples,
  c(
    "cases.submitter_id",
    "sample.submitter_id",
    "cases",
    "file_id",
    "file_name",
    "file_size",
    "sample_type",
    "is_ffpe"
  )
]

duplicate_exp_metadata
# ============================================================
# STEP 2C — Retrieve TCGA-BRCA biospecimen metadata
# ============================================================

biospecimen_raw <- pipeline_query(ctx, "biospecimen_raw", list(
  project = "TCGA-BRCA",
  type = "Biospecimen",
  save.csv = FALSE
), kind = "clinic")
dim(biospecimen_raw)

names(biospecimen_raw)

head(biospecimen_raw)
grep(
  "submitter|sample|aliquot|portion|analyte|vial|ffpe|preservation|tumor|weight|barcode",
  names(biospecimen_raw),
  value = TRUE,
  ignore.case = TRUE
)
# ============================================================
# STEP 2D — Match RNA-seq samples to GDC biospecimen metadata
# ============================================================

rna_sample_ids <- unique(
  exp_results$sample.submitter_id
)

length(rna_sample_ids)
# How many RNA-seq sample IDs are represented
# in the biospecimen metadata?
sum(
  rna_sample_ids %in%
    biospecimen_raw$submitter_id
)

# Which RNA-seq samples, if any, are missing
# from the biospecimen metadata?
rna_samples_missing_biospecimen <- setdiff(
  rna_sample_ids,
  biospecimen_raw$submitter_id
)

length(rna_samples_missing_biospecimen)

rna_samples_missing_biospecimen
rna_biospecimen <- biospecimen_raw[
  biospecimen_raw$submitter_id %in%
    rna_sample_ids,
]

dim(rna_biospecimen)

length(unique(rna_biospecimen$submitter_id))

anyDuplicated(rna_biospecimen$submitter_id)
table(
  rna_biospecimen$sample_type,
  useNA = "ifany"
)

table(
  rna_biospecimen$tissue_type,
  useNA = "ifany"
)

table(
  rna_biospecimen$tumor_descriptor,
  useNA = "ifany"
)

table(
  rna_biospecimen$preservation_method,
  useNA = "ifany"
)
# ============================================================
# Inspect biospecimen metadata for patients
# with multiple Primary Tumor RNA-seq samples
# ============================================================

duplicate_biospecimen <- rna_biospecimen[
  substr(rna_biospecimen$submitter_id, 1, 12) %in%
    patients_multiple_samples,
]

duplicate_biospecimen[
  order(duplicate_biospecimen$submitter_id),
  c(
    "submitter_id",
    "sample_type",
    "tumor_descriptor",
    "specimen_type",
    "tissue_type",
    "preservation_method",
    "initial_weight",
    "current_weight",
    "days_to_collection",
    "days_to_sample_procurement"
  )
]
# ============================================================
# STEP 2E — Inspect nested portion/analyte metadata
# ============================================================

class(biospecimen_raw$portions)

class(biospecimen_raw$portions[[1]])

str(
  biospecimen_raw$portions[[1]],
  max.level = 5
)
idx_A13E_01A <- which(
  biospecimen_raw$submitter_id ==
    "TCGA-A7-A13E-01A"
)

length(idx_A13E_01A)

str(
  biospecimen_raw$portions[[idx_A13E_01A]],
  max.level = 6
)
idx_A13E_01B <- which(
  biospecimen_raw$submitter_id ==
    "TCGA-A7-A13E-01B"
)

length(idx_A13E_01B)

str(
  biospecimen_raw$portions[[idx_A13E_01B]],
  max.level = 6
)
# ============================================================
# STEP 2F — Flatten RNA analyte metadata for multi-sample patients
# ============================================================

extract_rna_metadata <- function(sample_row) {
  
  sample_id <- rna_biospecimen$submitter_id[sample_row]
  preservation <- rna_biospecimen$preservation_method[sample_row]
  
  portions_df <- rna_biospecimen$portions[[sample_row]]
  
  # Return NULL if no portion information is available
  if (is.null(portions_df) || nrow(portions_df) == 0) {
    return(NULL)
  }
  
  output <- list()
  k <- 1
  
  for (p in seq_len(nrow(portions_df))) {
    
    analytes_df <- portions_df$analytes[[p]]
    
    if (is.null(analytes_df) || nrow(analytes_df) == 0) {
      next
    }
    
    # Keep RNA analytes only
    rna_idx <- which(analytes_df$analyte_type == "RNA")
    
    if (length(rna_idx) == 0) {
      next
    }
    
    for (a in rna_idx) {
      
      output[[k]] <- data.frame(
        sample_id = sample_id,
        preservation_method = preservation,
        portion_id = portions_df$submitter_id[p],
        portion_is_ffpe = portions_df$is_ffpe[p],
        analyte_id = analytes_df$submitter_id[a],
        analyte_type = analytes_df$analyte_type[a],
        rna_integrity_number =
          if ("rna_integrity_number" %in% names(analytes_df))
            analytes_df$rna_integrity_number[a]
        else NA,
        rrna_28s_18s_ratio =
          if ("ribosomal_rna_28s_18s_ratio" %in% names(analytes_df))
            analytes_df$ribosomal_rna_28s_18s_ratio[a]
        else NA,
        analyte_concentration =
          if ("concentration" %in% names(analytes_df))
            analytes_df$concentration[a]
        else NA,
        stringsAsFactors = FALSE
      )
      
      k <- k + 1
    }
  }
  
  if (length(output) == 0) {
    return(NULL)
  }
  
  do.call(rbind, output)
}
duplicate_sample_rows <- which(
  substr(rna_biospecimen$submitter_id, 1, 12) %in%
    patients_multiple_samples
)

duplicate_rna_quality <- do.call(
  rbind,
  lapply(
    duplicate_sample_rows,
    extract_rna_metadata
  )
)

duplicate_rna_quality[
  order(duplicate_rna_quality$sample_id),
]
# ============================================================
# STEP 2F validation — RNA quality by FFPE status
# ============================================================

# Archived exports can have an entirely missing optional QC metric.
# Skip only the descriptive summary, never impute or drop analysis samples.
diagnostic_aggregate <- function(formula, data, FUN) {
  observed <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  if (!nrow(observed)) {
    message("Optional descriptive QC summary unavailable: ", deparse(formula))
    return(invisible(NULL))
  }
  stats::aggregate(formula, data = data, FUN = FUN)
}

diagnostic_aggregate(
  rna_integrity_number ~ portion_is_ffpe,
  data = duplicate_rna_quality,
  FUN = function(x) c(
    n = length(x),
    min = min(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
)

diagnostic_aggregate(
  rrna_28s_18s_ratio ~ portion_is_ffpe,
  data = duplicate_rna_quality,
  FUN = function(x) c(
    n = length(x),
    min = min(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
)
# ============================================================
# STEP 2G — Apply specimen-level selection rule
# Keep non-FFPE sample when a patient has multiple samples
# ============================================================

selected_sample_ids <- rna_sample_ids

for (patient in patients_multiple_samples) {
  
  patient_samples <- duplicate_rna_quality[
    substr(duplicate_rna_quality$sample_id, 1, 12) == patient,
  ]
  
  preferred_sample <- patient_samples$sample_id[
    patient_samples$portion_is_ffpe == FALSE
  ]
  
  selected_sample_ids <- setdiff(
    selected_sample_ids,
    patient_samples$sample_id
  )
  
  selected_sample_ids <- c(
    selected_sample_ids,
    preferred_sample
  )
}

selected_sample_ids <- unique(selected_sample_ids)

length(selected_sample_ids)
exp_after_specimen_selection <- exp_results[
  exp_results$sample.submitter_id %in%
    selected_sample_ids,
]

nrow(exp_after_specimen_selection)

length(
  unique(exp_after_specimen_selection$sample.submitter_id)
)

length(
  unique(exp_after_specimen_selection$cases.submitter_id)
)
files_per_selected_sample <- table(
  exp_after_specimen_selection$sample.submitter_id
)

table(files_per_selected_sample)

technical_duplicate_samples <- names(
  files_per_selected_sample[
    files_per_selected_sample > 1
  ]
)

length(technical_duplicate_samples)

technical_duplicate_samples
technical_duplicate_metadata <-
  exp_after_specimen_selection[
    exp_after_specimen_selection$sample.submitter_id %in%
      technical_duplicate_samples,
    c(
      "cases.submitter_id",
      "sample.submitter_id",
      "cases",
      "file_id",
      "file_name",
      "analysis_id",
      "analysis_submitter_id",
      "created_datetime",
      "updated_datetime"
    )
  ]

technical_duplicate_metadata[
  order(
    technical_duplicate_metadata$sample.submitter_id,
    technical_duplicate_metadata$cases
  ),
]
# ============================================================
# STEP 2H — Extract RNA aliquot metadata for technical duplicates
# ============================================================

extract_rna_aliquots <- function(sample_row) {
  
  sample_id <- rna_biospecimen$submitter_id[sample_row]
  portions_df <- rna_biospecimen$portions[[sample_row]]
  
  if (is.null(portions_df) || nrow(portions_df) == 0) {
    return(NULL)
  }
  
  output <- list()
  k <- 1
  
  for (p in seq_len(nrow(portions_df))) {
    
    analytes_df <- portions_df$analytes[[p]]
    
    if (is.null(analytes_df) || nrow(analytes_df) == 0) {
      next
    }
    
    rna_idx <- which(
      analytes_df$analyte_type == "RNA"
    )
    
    if (length(rna_idx) == 0) {
      next
    }
    
    for (a in rna_idx) {
      
      aliquots_df <- analytes_df$aliquots[[a]]
      
      if (is.null(aliquots_df) || nrow(aliquots_df) == 0) {
        next
      }
      
      for (q in seq_len(nrow(aliquots_df))) {
        
        output[[k]] <- data.frame(
          
          sample_id = sample_id,
          
          portion_id =
            portions_df$submitter_id[p],
          
          analyte_id =
            analytes_df$submitter_id[a],
          
          rna_integrity_number =
            if ("rna_integrity_number" %in% names(analytes_df))
              analytes_df$rna_integrity_number[a]
          else NA,
          
          aliquot_id =
            aliquots_df$submitter_id[q],
          
          aliquot_quantity =
            if ("aliquot_quantity" %in% names(aliquots_df))
              aliquots_df$aliquot_quantity[q]
          else NA,
          
          aliquot_volume =
            if ("aliquot_volume" %in% names(aliquots_df))
              aliquots_df$aliquot_volume[q]
          else NA,
          
          concentration =
            if ("concentration" %in% names(aliquots_df))
              aliquots_df$concentration[q]
          else NA,
          
          source_center =
            if ("source_center" %in% names(aliquots_df))
              aliquots_df$source_center[q]
          else NA,
          
          stringsAsFactors = FALSE
        )
        
        k <- k + 1
      }
    }
  }
  
  if (length(output) == 0) {
    return(NULL)
  }
  
  do.call(rbind, output)
}
technical_sample_rows <- which(
  rna_biospecimen$submitter_id %in%
    technical_duplicate_samples
)

technical_aliquot_metadata <- do.call(
  rbind,
  lapply(
    technical_sample_rows,
    extract_rna_aliquots
  )
)
star_aliquot_ids <- unique(
  technical_duplicate_metadata$cases
)

technical_star_aliquots <-
  technical_aliquot_metadata[
    technical_aliquot_metadata$aliquot_id %in%
      star_aliquot_ids,
  ]

technical_star_aliquots[
  order(
    technical_star_aliquots$sample_id,
    technical_star_aliquots$aliquot_id
  ),
]
nrow(technical_star_aliquots)
length(unique(
  technical_star_aliquots$aliquot_id
))
# ============================================================
# STEP 2I — Prepare the 10 aliquot-level STAR-count files
# for expression-level QC
# ============================================================

technical_duplicate_file_ids <-
  technical_duplicate_metadata$file_id

length(technical_duplicate_file_ids)

length(unique(technical_duplicate_file_ids))
class(query_exp)

names(query_exp)

str(
  query_exp,
  max.level = 2
)
# ============================================================
# STEP 2I — Build a query object for the 10 technical-multiplicity files
# ============================================================

query_technical <- query_exp

query_technical$results[[1]] <-
  query_exp$results[[1]][
    query_exp$results[[1]]$file_id %in%
      technical_duplicate_file_ids,
  ]

dim(query_technical$results[[1]])

length(unique(
  query_technical$results[[1]]$file_id
))
setequal(
  query_technical$results[[1]]$file_id,
  technical_duplicate_file_ids
)
technical_query_results <- getResults(
  query_technical
)

dim(technical_query_results)

technical_query_results[
  ,
  c(
    "cases.submitter_id",
    "sample.submitter_id",
    "cases",
    "file_id",
    "file_name"
  )
]
technical_qc_dir <- file.path(input_dir, "data_raw", "rna_seq_technical_qc")
technical_count_files <- pipeline_count_files(
  ctx, query_technical, technical_query_results, "rna_seq_technical_qc", 10L
)

length(technical_count_files)

basename(technical_count_files)
# ============================================================
# STEP 2J/K — One STAR-count reader shared by technical and final cohorts
# ============================================================

read_star_counts <- function(file) {
  read.delim(
    file,
    skip = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

example_star <- read_star_counts(technical_count_files[1])
dim(example_star)
names(example_star)
head(example_star, 10)
str(example_star)
# Raw header inspection is diagnostic; parsing always uses skip = 1 above.
readLines(technical_count_files[1], n = 15)
count.fields(technical_count_files[1], sep = "\t")[1:15]
file.info(technical_count_files[1])[, "size"]

calculate_star_qc <- function(file) {
  
  dat <- read_star_counts(file)
  
  summary_ids <- c(
    "N_unmapped",
    "N_multimapping",
    "N_noFeature",
    "N_ambiguous"
  )
  
  summary_rows <- dat[
    dat$gene_id %in% summary_ids,
  ]
  
  gene_rows <- dat[
    grepl("^ENSG", dat$gene_id),
  ]
  
  data.frame(
    file_name = basename(file),
    
    n_gene_rows = nrow(gene_rows),
    
    assigned_counts =
      sum(gene_rows$unstranded, na.rm = TRUE),
    
    detected_genes =
      sum(gene_rows$unstranded > 0, na.rm = TRUE),
    
    genes_count_ge_10 =
      sum(gene_rows$unstranded >= 10, na.rm = TRUE),
    
    zero_fraction =
      mean(gene_rows$unstranded == 0, na.rm = TRUE),
    
    N_unmapped =
      summary_rows$unstranded[
        summary_rows$gene_id == "N_unmapped"
      ],
    
    N_multimapping =
      summary_rows$unstranded[
        summary_rows$gene_id == "N_multimapping"
      ],
    
    N_noFeature =
      summary_rows$unstranded[
        summary_rows$gene_id == "N_noFeature"
      ],
    
    N_ambiguous =
      summary_rows$unstranded[
        summary_rows$gene_id == "N_ambiguous"
      ],
    
    stringsAsFactors = FALSE
  )
}
technical_expression_qc <- do.call(
  rbind,
  lapply(
    technical_count_files,
    calculate_star_qc
  )
)

technical_expression_qc
technical_expression_qc <- merge(
  technical_expression_qc,
  technical_query_results[
    ,
    c(
      "file_name",
      "cases.submitter_id",
      "sample.submitter_id",
      "cases",
      "file_id"
    )
  ],
  by = "file_name",
  all.x = TRUE
)

technical_expression_qc <- technical_expression_qc[
  order(
    technical_expression_qc$cases.submitter_id,
    technical_expression_qc$cases
  ),
]

technical_expression_qc
nrow(technical_expression_qc)

table(
  technical_expression_qc$cases.submitter_id
)
# ============================================================
# STEP 2L — CANONICAL technical count list and matrix
# ============================================================

technical_count_list <- lapply(
  technical_count_files,
  function(file) {
    
    dat <- read_star_counts(file)
    
    gene_rows <- dat[
      grepl("^ENSG", dat$gene_id),
      c("gene_id", "unstranded")
    ]
    
    gene_rows
  }
)

names(technical_count_list) <-
  basename(technical_count_files)
stopifnot(
  length(technical_count_list) > 0L,
  !anyDuplicated(names(technical_count_list)),
  all(vapply(technical_count_list, function(x) {
    is.numeric(x$unstranded) && !anyNA(x$unstranded) &&
      all(is.finite(x$unstranded)) && all(x$unstranded >= 0) &&
      all(x$unstranded == round(x$unstranded)) &&
      !anyNA(x$gene_id) && all(nzchar(x$gene_id)) && !anyDuplicated(x$gene_id)
  }, logical(1)))
)
all_identical_gene_ids <- all(vapply(
  technical_count_list,
  function(x) identical(x$gene_id, technical_count_list[[1]]$gene_id),
  logical(1)
))
stopifnot(all_identical_gene_ids, nrow(technical_count_list[[1]]) > 1L)
names(technical_count_list[[1]])

length(technical_count_list[[1]]$gene_id)

length(technical_count_list[[1]]$unstranded)

sapply(
  technical_count_list,
  nrow
)
technical_count_matrix <- do.call(
  cbind,
  lapply(
    technical_count_list,
    function(x) x$unstranded
  )
)
rownames(technical_count_matrix) <-
  technical_count_list[[1]]$gene_id

colnames(technical_count_matrix) <-
  names(technical_count_list)
dim(technical_count_matrix)

class(technical_count_matrix)

length(rownames(technical_count_matrix))
length(colnames(technical_count_matrix))

anyDuplicated(rownames(technical_count_matrix))
anyDuplicated(colnames(technical_count_matrix))

anyNA(technical_count_matrix)
head(
  technical_count_matrix[, 1:3]
)
# ============================================================
# STEP 2M — Map matrix columns to patient and aliquot metadata
# ============================================================

technical_file_map <- technical_query_results[
  ,
  c(
    "file_name",
    "cases.submitter_id",
    "cases"
  )
]

technical_file_map <- technical_file_map[
  match(
    colnames(technical_count_matrix),
    technical_file_map$file_name
  ),
]

technical_file_map
identical(
  technical_file_map$file_name,
  colnames(technical_count_matrix)
)

anyNA(technical_file_map$file_name)

table(
  technical_file_map$cases.submitter_id
)
# ============================================================
# STEP 2N — Pairwise expression similarity across all genes
# ============================================================

log_counts <- log2(
  technical_count_matrix + 1
)


technical_patients <- unique(
  technical_file_map$cases.submitter_id
)
stopifnot(
  identical(technical_file_map$file_name, colnames(technical_count_matrix)),
  !anyNA(technical_file_map),
  !anyDuplicated(technical_file_map$file_name),
  all(table(technical_file_map$cases.submitter_id) == 2L)
)

pairwise_similarity <- do.call(
  rbind,
  lapply(
    technical_patients,
    function(patient) {
      
      idx <- which(
        technical_file_map$cases.submitter_id == patient
      )
      
      data.frame(
        patient_id = patient,
        
        aliquot_1 =
          technical_file_map$cases[idx[1]],
        
        aliquot_2 =
          technical_file_map$cases[idx[2]],
        
        pearson_log_counts =
          cor(
            log_counts[, idx[1]],
            log_counts[, idx[2]],
            method = "pearson",
            use = "complete.obs"
          ),
        
        spearman_counts =
          cor(
            technical_count_matrix[, idx[1]],
            technical_count_matrix[, idx[2]],
            method = "spearman",
            use = "complete.obs"
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)

pairwise_similarity
# ============================================================
# STEP 2O — Similarity among genes expressed in at least
# one member of each aliquot pair
# ============================================================

pairwise_similarity_expressed <- do.call(
  rbind,
  lapply(
    technical_patients,
    function(patient) {
      
      idx <- which(
        technical_file_map$cases.submitter_id == patient
      )
      
      c1 <- technical_count_matrix[, idx[1]]
      c2 <- technical_count_matrix[, idx[2]]
      
      keep <- c1 > 0 | c2 > 0
      
      data.frame(
        patient_id = patient,
        
        n_genes_compared =
          sum(keep),
        
        pearson_expressed =
          cor(
            log2(c1[keep] + 1),
            log2(c2[keep] + 1),
            method = "pearson"
          ),
        
        spearman_expressed =
          cor(
            c1[keep],
            c2[keep],
            method = "spearman"
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)

pairwise_similarity_expressed
pairwise_depth_comparison <- do.call(
  rbind,
  lapply(
    technical_patients,
    function(patient) {
      
      idx <- which(
        technical_file_map$cases.submitter_id == patient
      )
      
      lib1 <- sum(
        technical_count_matrix[, idx[1]]
      )
      
      lib2 <- sum(
        technical_count_matrix[, idx[2]]
      )
      
      data.frame(
        patient_id = patient,
        
        aliquot_1 =
          technical_file_map$cases[idx[1]],
        
        aliquot_2 =
          technical_file_map$cases[idx[2]],
        
        assigned_counts_1 = lib1,
        assigned_counts_2 = lib2,
        
        depth_ratio =
          max(lib1, lib2) / min(lib1, lib2),
        
        stringsAsFactors = FALSE
      )
    }
  )
)

pairwise_depth_comparison
# ============================================================
# STEP 2P — Library-size adjusted aliquot similarity
# Diagnostic QC only
# ============================================================

technical_cpm <- sweep(
  technical_count_matrix,
  2,
  colSums(technical_count_matrix),
  FUN = "/"
) * 1e6

log2_cpm <- log2(
  technical_cpm + 1
)
dim(technical_cpm)

summary(
  colSums(technical_cpm)
)
pairwise_similarity_cpm <- do.call(
  rbind,
  lapply(
    technical_patients,
    function(patient) {
      
      idx <- which(
        technical_file_map$cases.submitter_id == patient
      )
      
      cpm1 <- technical_cpm[, idx[1]]
      cpm2 <- technical_cpm[, idx[2]]
      
      keep <- cpm1 >= 1 | cpm2 >= 1
      
      data.frame(
        patient_id = patient,
        
        aliquot_1 =
          technical_file_map$cases[idx[1]],
        
        aliquot_2 =
          technical_file_map$cases[idx[2]],
        
        n_genes_compared =
          sum(keep),
        
        pearson_log2_cpm =
          cor(
            log2_cpm[keep, idx[1]],
            log2_cpm[keep, idx[2]],
            method = "pearson"
          ),
        
        spearman_cpm =
          cor(
            cpm1[keep],
            cpm2[keep],
            method = "spearman"
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)

pairwise_similarity_cpm
technical_pair_summary <- merge(
  pairwise_similarity,
  pairwise_similarity_expressed,
  by = "patient_id"
)

technical_pair_summary <- merge(
  technical_pair_summary,
  pairwise_depth_comparison[
    ,
    c(
      "patient_id",
      "depth_ratio"
    )
  ],
  by = "patient_id"
)

technical_pair_summary <- merge(
  technical_pair_summary,
  pairwise_similarity_cpm[
    ,
    c(
      "patient_id",
      "pearson_log2_cpm",
      "spearman_cpm"
    )
  ],
  by = "patient_id"
)

technical_pair_summary
technical_pair_summary[
  ,
  c(
    "patient_id",
    "pearson_log_counts",
    "pearson_expressed",
    "pearson_log2_cpm",
    "spearman_counts",
    "spearman_expressed",
    "spearman_cpm",
    "depth_ratio"
  )
]
summary(colSums(technical_cpm))

pairwise_similarity_cpm

technical_pair_summary[
  ,
  c(
    "patient_id",
    "pearson_log_counts",
    "pearson_expressed",
    "pearson_log2_cpm",
    "spearman_counts",
    "spearman_expressed",
    "spearman_cpm",
    "depth_ratio"
  )
]
# ============================================================
# STEP 2Q — Gene-assignment fraction for paired aliquots
# ============================================================

technical_expression_qc$total_star_accounted <-
  technical_expression_qc$assigned_counts +
  technical_expression_qc$N_unmapped +
  technical_expression_qc$N_multimapping +
  technical_expression_qc$N_noFeature +
  technical_expression_qc$N_ambiguous

technical_expression_qc$gene_assignment_fraction <-
  technical_expression_qc$assigned_counts /
  technical_expression_qc$total_star_accounted

technical_expression_qc$gene_assignment_percent <-
  technical_expression_qc$gene_assignment_fraction * 100
technical_assignment_qc <- technical_expression_qc[
  ,
  c(
    "cases.submitter_id",
    "cases",
    "assigned_counts",
    "detected_genes",
    "zero_fraction",
    "N_unmapped",
    "N_multimapping",
    "N_noFeature",
    "N_ambiguous",
    "gene_assignment_percent"
  )
]

technical_assignment_qc <- technical_assignment_qc[
  order(
    technical_assignment_qc$cases.submitter_id,
    technical_assignment_qc$cases
  ),
]

technical_assignment_qc
# ============================================================
# STEP 2R — Final aliquot-level selection
# ============================================================

selected_technical_aliquots <- c(
  "TCGA-A7-A0DB-01A-11R-A00Z-07",
  "TCGA-A7-A13D-01A-13R-A12P-07",
  "TCGA-A7-A13E-01A-11R-A12P-07",
  "TCGA-A7-A26E-01A-11R-A169-07",
  "TCGA-A7-A26J-01A-11R-A169-07"
)

selected_technical_aliquots
technical_duplicate_patients <- c(
  "TCGA-A7-A0DB",
  "TCGA-A7-A13D",
  "TCGA-A7-A13E",
  "TCGA-A7-A26E",
  "TCGA-A7-A26J"
)
keep_final_file <-
  !exp_after_specimen_selection$cases.submitter_id %in%
  technical_duplicate_patients |
  exp_after_specimen_selection$cases %in%
  selected_technical_aliquots

exp_final_selected <- exp_after_specimen_selection[
  keep_final_file,
]
nrow(exp_final_selected)

length(
  unique(exp_final_selected$file_id)
)

length(
  unique(exp_final_selected$cases)
)

length(
  unique(exp_final_selected$sample.submitter_id)
)

length(
  unique(exp_final_selected$cases.submitter_id)
)
table(
  table(
    exp_final_selected$cases.submitter_id
  )
)
exp_final_selected[
  exp_final_selected$cases.submitter_id %in%
    technical_duplicate_patients,
  c(
    "cases.submitter_id",
    "sample.submitter_id",
    "cases",
    "file_id",
    "file_name"
  )
]
rna_patient_ids <- unique(
  exp_final_selected$cases.submitter_id
)

clinical_patient_ids <- unique(
  clinical_clean$patient_id
)

length(rna_patient_ids)
length(clinical_patient_ids)

length(
  intersect(
    rna_patient_ids,
    clinical_patient_ids
  )
)

setdiff(
  rna_patient_ids,
  clinical_patient_ids
)

setdiff(
  clinical_patient_ids,
  rna_patient_ids
)
# ============================================================
pipeline_assert_ids(exp_final_selected$cases.submitter_id, "selected RNA patients")
pipeline_assert_ids(exp_final_selected$file_id, "selected RNA files")
stopifnot(all(exp_final_selected$sample_type == "Primary Tumor"))
# STEP 2S — Build final matched molecular-clinical cohort
# ============================================================

rna_metadata_final <- exp_final_selected[
  exp_final_selected$cases.submitter_id %in%
    clinical_patient_ids,
  c(
    "cases.submitter_id",
    "sample.submitter_id",
    "cases",
    "file_id",
    "file_name"
  )
]
names(rna_metadata_final)[
  names(rna_metadata_final) == "cases.submitter_id"
] <- "patient_id"
matched_cohort <- merge(
  clinical_clean,
  rna_metadata_final,
  by = "patient_id",
  all = FALSE
)
dim(matched_cohort)

length(unique(matched_cohort$patient_id))

anyDuplicated(matched_cohort$patient_id)

anyNA(matched_cohort$file_id)

table(matched_cohort$age_available)

table(matched_cohort$stage_available)

table(matched_cohort$age_stage_available)
data.frame(
  cohort_stage = c(
    "RNA-seq patients",
    "Clinical curated patients",
    "Matched molecular-clinical patients",
    "Age-analysis eligible",
    "Stage-analysis eligible",
    "Age+Stage available"
  ),
  
  n = c(
    length(rna_patient_ids),
    length(clinical_patient_ids),
    nrow(matched_cohort),
    sum(matched_cohort$age_available),
    sum(matched_cohort$stage_available),
    sum(matched_cohort$age_stage_available)
  )
)
pipeline_assert_ids(matched_cohort$patient_id, "matched patients")
pipeline_assert_ids(matched_cohort$file_id, "matched files")
stopifnot(all(matched_cohort$patient_id %in% clinical_clean$patient_id))
write.csv(
  matched_cohort,
  file.path(
    output_dir,
    "data_processed",
    "TCGA_BRCA_matched_molecular_clinical_cohort.csv"
  ),
  row.names = FALSE
)

saveRDS(
  matched_cohort,
  file.path(
    output_dir,
    "data_processed",
    "TCGA_BRCA_matched_molecular_clinical_cohort.rds"
  )
)
matched_cohort_check <- readRDS(
  file.path(
    output_dir,
    "data_processed",
    "TCGA_BRCA_matched_molecular_clinical_cohort.rds"
  )
)

stopifnot(identical(matched_cohort, matched_cohort_check))

# ============================================================
pipeline_event(ctx, "counts_acquisition")
# STEP 3A — Define the final RNA-seq file set
# ============================================================

final_rna_file_ids <- matched_cohort$file_id

length(final_rna_file_ids)

length(unique(final_rna_file_ids))

anyNA(final_rna_file_ids)
all(
  final_rna_file_ids %in%
    query_exp$results[[1]]$file_id
)
sum(
  query_exp$results[[1]]$file_id %in%
    final_rna_file_ids
)
# ============================================================
# STEP 3A — Construct final GDC query
# ============================================================

query_final_rna <- query_exp

query_final_rna$results[[1]] <-
  query_exp$results[[1]][
    query_exp$results[[1]]$file_id %in%
      final_rna_file_ids,
  ]
final_query_results <- getResults(
  query_final_rna
)

dim(final_query_results)
nrow(final_query_results)

length(
  unique(final_query_results$file_id)
)

length(
  unique(final_query_results$cases.submitter_id)
)

length(
  unique(final_query_results$sample.submitter_id)
)

length(
  unique(final_query_results$cases)
)
setequal(
  final_query_results$file_id,
  matched_cohort$file_id
)
setdiff(
  matched_cohort$file_id,
  final_query_results$file_id
)

setdiff(
  final_query_results$file_id,
  matched_cohort$file_id
)
setequal(
  final_query_results$cases.submitter_id,
  matched_cohort$patient_id
)
table(
  table(
    final_query_results$cases.submitter_id
  )
)
# ============================================================
# STEP 3A — Create reproducible RNA-seq manifest
# ============================================================

final_rna_manifest <- final_query_results[
  ,
  c(
    "cases.submitter_id",
    "sample.submitter_id",
    "cases",
    "file_id",
    "file_name",
    "md5sum",
    "file_size",
    "data_type",
    "analysis_workflow_type"
  )
]
names(final_rna_manifest)[
  names(final_rna_manifest) == "cases.submitter_id"
] <- "patient_id"

names(final_rna_manifest)[
  names(final_rna_manifest) == "cases"
] <- "aliquot_id"
final_rna_manifest <- final_rna_manifest[
  match(
    matched_cohort$patient_id,
    final_rna_manifest$patient_id
  ),
]
identical(
  final_rna_manifest$patient_id,
  matched_cohort$patient_id
)

nrow(final_rna_manifest)

anyDuplicated(final_rna_manifest$patient_id)

anyDuplicated(final_rna_manifest$file_id)

anyNA(final_rna_manifest$file_id)
pipeline_assert_ids(final_rna_manifest$patient_id, "manifest patients")
pipeline_assert_ids(final_rna_manifest$file_id, "manifest files")
pipeline_assert_ids(final_rna_manifest$file_name, "manifest filenames")
stopifnot(identical(final_rna_manifest$patient_id, matched_cohort$patient_id),
          setequal(final_rna_manifest$file_id, matched_cohort$file_id),
          !anyNA(final_rna_manifest$md5sum))
write.csv(
  final_rna_manifest,
  file.path(
    output_dir,
    "data_processed",
    "TCGA_BRCA_final_RNAseq_manifest.csv"
  ),
  row.names = FALSE
)

saveRDS(
  final_rna_manifest,
  file.path(
    output_dir,
    "data_processed",
    "TCGA_BRCA_final_RNAseq_manifest.rds"
  )
)

length(final_rna_file_ids)

length(unique(final_rna_file_ids))

anyNA(final_rna_file_ids)

all(
  final_rna_file_ids %in%
    query_exp$results[[1]]$file_id
)

dim(final_query_results)

setequal(
  final_query_results$file_id,
  matched_cohort$file_id
)

setequal(
  final_query_results$cases.submitter_id,
  matched_cohort$patient_id
)

table(
  table(
    final_query_results$cases.submitter_id
  )
)

identical(
  final_rna_manifest$patient_id,
  matched_cohort$patient_id
)

# ============================================================
# STEP 3B — Download size audit
# ============================================================

total_download_bytes <- sum(
  final_query_results$file_size,
  na.rm = TRUE
)

total_download_gb <- total_download_bytes / 1024^3

total_download_gb
summary(
  final_query_results$file_size / 1024^2
)
data.frame(
  n_files = nrow(final_query_results),
  total_GB = total_download_gb,
  mean_MB = mean(final_query_results$file_size) / 1024^2,
  median_MB = median(final_query_results$file_size) / 1024^2,
  max_MB = max(final_query_results$file_size) / 1024^2
)
# ============================================================
# STEP 3C — Download final RNA-seq cohort
# ============================================================

final_rna_dir <- file.path(input_dir, "data_raw", "rna_seq_final")
downloaded_star_files <- pipeline_count_files(
  ctx, query_final_rna, final_query_results, "rna_seq_final", 50L
)

length(downloaded_star_files)
downloaded_file_names <- basename(
  downloaded_star_files
)

length(unique(downloaded_file_names))

setequal(
  downloaded_file_names,
  final_rna_manifest$file_name
)
setdiff(
  final_rna_manifest$file_name,
  downloaded_file_names
)

setdiff(
  downloaded_file_names,
  final_rna_manifest$file_name
)
anyDuplicated(
  downloaded_file_names
)
length(downloaded_star_files)

length(unique(downloaded_file_names))

setequal(
  downloaded_file_names,
  final_rna_manifest$file_name
)

setdiff(
  final_rna_manifest$file_name,
  downloaded_file_names
)

setdiff(
  downloaded_file_names,
  final_rna_manifest$file_name
)

anyDuplicated(downloaded_file_names)

# ============================================================
# STEP 3E — Validate downloaded files using MD5 checksums
# ============================================================

downloaded_md5 <- tools::md5sum(
  downloaded_star_files
)
length(downloaded_md5)

head(downloaded_md5)
md5_validation <- data.frame(
  file_name = basename(names(downloaded_md5)),
  md5_observed = unname(downloaded_md5),
  stringsAsFactors = FALSE
)
md5_validation$md5_expected <-
  final_rna_manifest$md5sum[
    match(
      md5_validation$file_name,
      final_rna_manifest$file_name
    )
  ]
md5_validation$md5_match <-
  tolower(md5_validation$md5_observed) ==
  tolower(md5_validation$md5_expected)
nrow(md5_validation)

anyNA(md5_validation$md5_expected)

table(md5_validation$md5_match)
md5_validation[
  !md5_validation$md5_match |
    is.na(md5_validation$md5_match),
]
data.frame(
  files_checked = nrow(md5_validation),
  files_passed = sum(md5_validation$md5_match, na.rm = TRUE),
  files_failed = sum(!md5_validation$md5_match, na.rm = TRUE),
  missing_expected_md5 = sum(is.na(md5_validation$md5_expected))
)
# ============================================================
stopifnot(nrow(md5_validation) == nrow(final_rna_manifest),
          !anyNA(md5_validation$md5_match), all(md5_validation$md5_match))
# STEP 3F — Parse final STAR-count cohort
# ============================================================

star_structure_qc <- lapply(
  downloaded_star_files,
  function(file) {
    
    dat <- read_star_counts(file)
    
    data.frame(
      file_name = basename(file),
      n_rows = nrow(dat),
      n_cols = ncol(dat),
      has_gene_id = "gene_id" %in% names(dat),
      has_unstranded = "unstranded" %in% names(dat),
      stringsAsFactors = FALSE
    )
  }
)

star_structure_qc <- do.call(
  rbind,
  star_structure_qc
)
summary(star_structure_qc$n_rows)

summary(star_structure_qc$n_cols)

table(star_structure_qc$has_gene_id)

table(star_structure_qc$has_unstranded)
table(
  star_structure_qc$n_rows
)

table(
  star_structure_qc$n_cols
)

final_count_list <- lapply(
  downloaded_star_files,
  function(file) {
    
    dat <- read_star_counts(file)
    
    gene_rows <- dat[
      grepl("^ENSG", dat$gene_id),
      c("gene_id", "unstranded")
    ]
    
    gene_rows
  }
)

names(final_count_list) <- basename(
  downloaded_star_files
)
length(final_count_list)

summary(
  vapply(
    final_count_list,
    nrow,
    integer(1)
  )
)
reference_gene_ids <- final_count_list[[1]]$gene_id

gene_id_consistency <- vapply(
  final_count_list,
  function(x) {
    identical(
      x$gene_id,
      reference_gene_ids
    )
  },
  logical(1)
)

table(gene_id_consistency)
length(reference_gene_ids)

anyDuplicated(reference_gene_ids)


anyNA(reference_gene_ids)
summary(star_structure_qc$n_rows)
summary(star_structure_qc$n_cols)
table(star_structure_qc$has_gene_id)
table(star_structure_qc$has_unstranded)

length(final_count_list)
summary(vapply(final_count_list, nrow, integer(1)))

table(gene_id_consistency)

length(reference_gene_ids)
anyDuplicated(reference_gene_ids)
anyNA(reference_gene_ids)
# ============================================================
stopifnot(all(gene_id_consistency), !anyNA(reference_gene_ids),
          !anyDuplicated(reference_gene_ids), all(star_structure_qc$has_gene_id),
          all(star_structure_qc$has_unstranded))
# STEP 3H — Construct canonical raw count matrix
# ============================================================

raw_count_matrix <- do.call(
  cbind,
  lapply(
    final_count_list,
    function(x) x$unstranded
  )
)

rownames(raw_count_matrix) <- reference_gene_ids

colnames(raw_count_matrix) <- names(
  final_count_list
)
dim(raw_count_matrix)

class(raw_count_matrix)

length(rownames(raw_count_matrix))

length(colnames(raw_count_matrix))
anyNA(raw_count_matrix)

anyDuplicated(
  rownames(raw_count_matrix)
)

anyDuplicated(
  colnames(raw_count_matrix)
)

all(
  raw_count_matrix >= 0
)
all(
  raw_count_matrix ==
    floor(raw_count_matrix)
)
matrix_file_names <- colnames(
  raw_count_matrix
)

matrix_patient_ids <-
  final_rna_manifest$patient_id[
    match(
      matrix_file_names,
      final_rna_manifest$file_name
    )
  ]
length(matrix_patient_ids)

anyNA(matrix_patient_ids)

anyDuplicated(matrix_patient_ids)
colnames(raw_count_matrix) <-
  matrix_patient_ids
length(unique(
  colnames(raw_count_matrix)
))
setequal(
  colnames(raw_count_matrix),
  matched_cohort$patient_id
)
raw_count_matrix <- raw_count_matrix[
  ,
  matched_cohort$patient_id,
  drop = FALSE
]
identical(
  colnames(raw_count_matrix),
  matched_cohort$patient_id
)
dim(raw_count_matrix)

typeof(raw_count_matrix)

range(raw_count_matrix)

summary(
  colSums(raw_count_matrix)
)
summary(
  rowSums(raw_count_matrix)
)
sum(
  rowSums(raw_count_matrix) == 0
)
dim(raw_count_matrix)

class(raw_count_matrix)

anyNA(raw_count_matrix)

anyDuplicated(rownames(raw_count_matrix))

anyDuplicated(colnames(raw_count_matrix))

all(raw_count_matrix >= 0)

all(
  raw_count_matrix ==
    floor(raw_count_matrix)
)

length(matrix_patient_ids)

anyNA(matrix_patient_ids)

anyDuplicated(matrix_patient_ids)

setequal(
  colnames(raw_count_matrix),
  matched_cohort$patient_id
)

identical(
  colnames(raw_count_matrix),
  matched_cohort$patient_id
)

summary(
  colSums(raw_count_matrix)
)

sum(
  rowSums(raw_count_matrix) == 0
)
pipeline_assert_counts(raw_count_matrix)
stopifnot(identical(colnames(raw_count_matrix), matched_cohort$patient_id))
saveRDS(
  raw_count_matrix,
  file.path(
    output_dir,
    "data_processed",
    "TCGA_BRCA_raw_count_matrix.rds"
  )
)

raw_count_matrix_check <- readRDS(
  file.path(
    output_dir,
    "data_processed",
    "TCGA_BRCA_raw_count_matrix.rds"
  )
)

stopifnot(identical(raw_count_matrix, raw_count_matrix_check))
# ============================================================
pipeline_event(ctx, "raw_count_QC")
# STEP 4A — Sample-level raw-count QC
# ============================================================

sample_qc <- data.frame(
  patient_id = colnames(raw_count_matrix),
  
  assigned_gene_counts =
    colSums(raw_count_matrix),
  
  detected_genes =
    colSums(raw_count_matrix > 0),
  
  genes_count_ge_10 =
    colSums(raw_count_matrix >= 10),
  
  zero_fraction =
    colMeans(raw_count_matrix == 0),
  
  stringsAsFactors = FALSE
)

sample_qc <- merge(
  sample_qc,
  matched_cohort[
    ,
    c(
      "patient_id",
      "age_years",
      "stage_group",
      "sex"
    )
  ],
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)
sample_qc <- sample_qc[
  match(
    colnames(raw_count_matrix),
    sample_qc$patient_id
  ),
]
identical(
  sample_qc$patient_id,
  colnames(raw_count_matrix)
)

anyDuplicated(sample_qc$patient_id)
summary(
  sample_qc$assigned_gene_counts
)

summary(
  sample_qc$detected_genes
)

summary(
  sample_qc$genes_count_ge_10
)

summary(
  sample_qc$zero_fraction
)
cor(
  sample_qc$assigned_gene_counts,
  sample_qc$detected_genes,
  method = "spearman"
)
sample_qc[
  order(sample_qc$assigned_gene_counts),
  c(
    "patient_id",
    "assigned_gene_counts",
    "detected_genes",
    "zero_fraction"
  )
][1:10,]
sample_qc[
  order(sample_qc$assigned_gene_counts, decreasing = TRUE),
  c(
    "patient_id",
    "assigned_gene_counts",
    "detected_genes",
    "zero_fraction"
  )
][1:10,]
sample_qc[
  order(sample_qc$detected_genes),
  c(
    "patient_id",
    "assigned_gene_counts",
    "detected_genes",
    "zero_fraction"
  )
][1:10,]
sample_qc[
  order(sample_qc$zero_fraction, decreasing = TRUE),
  c(
    "patient_id",
    "assigned_gene_counts",
    "detected_genes",
    "zero_fraction"
  )
][1:10,]

# ============================================================
# STEP 4B — Prepare sample-level QC diagnostics
# ============================================================

sample_qc$log10_assigned_counts <-
  log10(sample_qc$assigned_gene_counts)
library(ggplot2)

p_library <- ggplot(
  sample_qc,
  aes(x = assigned_gene_counts / 1e6)
) +
  geom_histogram(
    bins = 40,
    color = "white"
  ) +
  geom_vline(
    xintercept =
      median(sample_qc$assigned_gene_counts) / 1e6,
    linetype = "dashed"
  ) +
  labs(
    title = "Distribution of assigned gene counts",
    x = "Assigned gene counts (millions)",
    y = "Number of samples"
  ) +
  theme_minimal()

p_library
p_detected <- ggplot(
  sample_qc,
  aes(x = detected_genes)
) +
  geom_histogram(
    bins = 40,
    color = "white"
  ) +
  geom_vline(
    xintercept =
      median(sample_qc$detected_genes),
    linetype = "dashed"
  ) +
  labs(
    title = "Distribution of detected genes",
    x = "Genes with raw count > 0",
    y = "Number of samples"
  ) +
  theme_minimal()

p_detected
p_depth_detected <- ggplot(
  sample_qc,
  aes(
    x = assigned_gene_counts / 1e6,
    y = detected_genes
  )
) +
  geom_point(
    alpha = 0.55
  ) +
  geom_smooth(
    method = "loess",
    se = FALSE
  ) +
  labs(
    title = "Library size versus detected genes",
    subtitle = "Spearman rho = 0.608",
    x = "Assigned gene counts (millions)",
    y = "Detected genes"
  ) +
  theme_minimal()

p_depth_detected
p_zero <- ggplot(
  sample_qc,
  aes(x = zero_fraction)
) +
  geom_histogram(
    bins = 40,
    color = "white"
  ) +
  geom_vline(
    xintercept =
      median(sample_qc$zero_fraction),
    linetype = "dashed"
  ) +
  labs(
    title = "Distribution of zero-count fraction",
    x = "Fraction of genes with zero counts",
    y = "Number of samples"
  ) +
  theme_minimal()

p_zero
library_iqr <- IQR(
  sample_qc$log10_assigned_counts
)

library_q1 <- quantile(
  sample_qc$log10_assigned_counts,
  0.25
)

library_q3 <- quantile(
  sample_qc$log10_assigned_counts,
  0.75
)

library_lower <-
  library_q1 - 1.5 * library_iqr

library_upper <-
  library_q3 + 1.5 * library_iqr
detected_iqr <- IQR(
  sample_qc$detected_genes
)

detected_q1 <- quantile(
  sample_qc$detected_genes,
  0.25
)

detected_q3 <- quantile(
  sample_qc$detected_genes,
  0.75
)

detected_lower <-
  detected_q1 - 1.5 * detected_iqr

detected_upper <-
  detected_q3 + 1.5 * detected_iqr
sample_qc$library_flag <-
  sample_qc$log10_assigned_counts < library_lower |
  sample_qc$log10_assigned_counts > library_upper

sample_qc$detected_flag <-
  sample_qc$detected_genes < detected_lower |
  sample_qc$detected_genes > detected_upper
sample_qc$raw_qc_flag <-
  sample_qc$library_flag |
  sample_qc$detected_flag
table(sample_qc$library_flag)

table(sample_qc$detected_flag)

table(sample_qc$raw_qc_flag)
raw_qc_flagged_samples <- sample_qc[
  sample_qc$raw_qc_flag,
  c(
    "patient_id",
    "assigned_gene_counts",
    "detected_genes",
    "zero_fraction",
    "library_flag",
    "detected_flag",
    "age_years",
    "stage_group",
    "sex"
  )
]

raw_qc_flagged_samples <- raw_qc_flagged_samples[
  order(
    raw_qc_flagged_samples$detected_genes
  ),
]

raw_qc_flagged_samples
sample_qc[
  which.max(sample_qc$detected_genes),
  c(
    "patient_id",
    "assigned_gene_counts",
    "detected_genes",
    "genes_count_ge_10",
    "zero_fraction",
    "age_years",
    "stage_group",
    "sex"
  )
]

# ============================================================
# STEP 4C — STAR-level technical QC
# ============================================================

star_file_map <- data.frame(
  file_name = basename(downloaded_star_files),
  local_path = downloaded_star_files,
  stringsAsFactors = FALSE
)

star_file_map$patient_id <-
  final_rna_manifest$patient_id[
    match(
      star_file_map$file_name,
      final_rna_manifest$file_name
    )
  ]

stopifnot(
  !anyNA(star_file_map$patient_id),
  anyDuplicated(star_file_map$patient_id) == 0
)
test_star <- read_star_counts(
  star_file_map$local_path[1]
)

test_star[
  1:4,
  c("gene_id", "unstranded")
]
star_summary_list <- lapply(
  seq_len(nrow(star_file_map)),
  function(i) {
    
    dat <- read_star_counts(
      star_file_map$local_path[i]
    )
    
    summary_rows <- dat[
      1:4,
      c("gene_id", "unstranded")
    ]
    
    assigned_counts <- sum(
      dat$unstranded[
        grepl("^ENSG", dat$gene_id)
      ]
    )
    
    values <- setNames(
      summary_rows$unstranded,
      summary_rows$gene_id
    )
    
    data.frame(
      patient_id = star_file_map$patient_id[i],
      assigned_counts = assigned_counts,
      N_unmapped = unname(values["N_unmapped"]),
      N_multimapping = unname(values["N_multimapping"]),
      N_noFeature = unname(values["N_noFeature"]),
      N_ambiguous = unname(values["N_ambiguous"]),
      stringsAsFactors = FALSE
    )
  }
)

star_qc <- do.call(
  rbind,
  star_summary_list
)
assigned_match <- identical(
  star_qc$assigned_counts[
    match(
      sample_qc$patient_id,
      star_qc$patient_id
    )
  ],
  sample_qc$assigned_gene_counts
)

assigned_match
star_qc$total_accounted_counts <-
  star_qc$assigned_counts +
  star_qc$N_unmapped +
  star_qc$N_multimapping +
  star_qc$N_noFeature +
  star_qc$N_ambiguous

star_qc$gene_assignment_fraction <-
  star_qc$assigned_counts /
  star_qc$total_accounted_counts

star_qc$unmapped_fraction <-
  star_qc$N_unmapped /
  star_qc$total_accounted_counts

star_qc$multimapping_fraction <-
  star_qc$N_multimapping /
  star_qc$total_accounted_counts

star_qc$noFeature_fraction <-
  star_qc$N_noFeature /
  star_qc$total_accounted_counts

star_qc$ambiguous_fraction <-
  star_qc$N_ambiguous /
  star_qc$total_accounted_counts
sample_qc_star <- merge(
  sample_qc,
  star_qc,
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)

sample_qc_star <- sample_qc_star[
  match(
    sample_qc$patient_id,
    sample_qc_star$patient_id
  ),
]

stopifnot(
  identical(
    sample_qc_star$patient_id,
    sample_qc$patient_id
  )
)
summary(
  sample_qc_star$gene_assignment_fraction
)

summary(
  sample_qc_star$unmapped_fraction
)

summary(
  sample_qc_star$multimapping_fraction
)

summary(
  sample_qc_star$noFeature_fraction
)

summary(
  sample_qc_star$ambiguous_fraction
)
cor(
  sample_qc_star$gene_assignment_fraction,
  sample_qc_star$detected_genes,
  method = "spearman"
)

cor(
  sample_qc_star$gene_assignment_fraction,
  sample_qc_star$zero_fraction,
  method = "spearman"
)
flagged_star_qc <- sample_qc_star[
  sample_qc_star$raw_qc_flag,
  c(
    "patient_id",
    "assigned_gene_counts",
    "detected_genes",
    "zero_fraction",
    "gene_assignment_fraction",
    "unmapped_fraction",
    "multimapping_fraction",
    "noFeature_fraction",
    "ambiguous_fraction",
    "library_flag",
    "detected_flag"
  )
]

flagged_star_qc <- flagged_star_qc[
  order(
    flagged_star_qc$detected_genes
  ),
]

head(
  flagged_star_qc,
  15
)

tail(
  flagged_star_qc,
  15
)
sample_qc_star[
  order(
    sample_qc_star$detected_genes,
    decreasing = TRUE
  ),
  c(
    "patient_id",
    "assigned_gene_counts",
    "detected_genes",
    "genes_count_ge_10",
    "zero_fraction",
    "gene_assignment_fraction",
    "unmapped_fraction",
    "multimapping_fraction",
    "noFeature_fraction",
    "ambiguous_fraction"
  )
][1:15, ]
previous_multi_vial_patients <- c(
  "TCGA-A7-A0DB",
  "TCGA-A7-A0DC",
  "TCGA-A7-A13D",
  "TCGA-A7-A13E",
  "TCGA-A7-A13G",
  "TCGA-A7-A26E",
  "TCGA-A7-A26F",
  "TCGA-A7-A26I",
  "TCGA-A7-A26J",
  "TCGA-AC-A2QH",
  "TCGA-AC-A3OD"
)

sample_qc_star[
  sample_qc_star$patient_id %in%
    previous_multi_vial_patients,
  c(
    "patient_id",
    "assigned_gene_counts",
    "detected_genes",
    "zero_fraction",
    "gene_assignment_fraction",
    "raw_qc_flag"
  )
]
# ============================================================
# STEP 4C.1 — Diagnose assigned-count validation mismatch
# ============================================================

star_assigned_aligned <-
  star_qc$assigned_counts[
    match(
      sample_qc$patient_id,
      star_qc$patient_id
    )
  ]

length(star_assigned_aligned)
length(sample_qc$assigned_gene_counts)

typeof(star_assigned_aligned)
typeof(sample_qc$assigned_gene_counts)

all.equal(
  star_assigned_aligned,
  sample_qc$assigned_gene_counts
)

all(
  star_assigned_aligned ==
    sample_qc$assigned_gene_counts
)

max(
  abs(
    star_assigned_aligned -
      sample_qc$assigned_gene_counts
  )
)

sum(
  star_assigned_aligned !=
    sample_qc$assigned_gene_counts
)
assigned_value_match <- all(
  star_assigned_aligned ==
    sample_qc$assigned_gene_counts
)

assigned_value_match
# ============================================================
# STEP 4D — STAR technical outlier screening
# ============================================================

star_metrics <- c(
  "gene_assignment_fraction",
  "unmapped_fraction",
  "multimapping_fraction",
  "noFeature_fraction",
  "ambiguous_fraction"
)

star_iqr_thresholds <- do.call(
  rbind,
  lapply(
    star_metrics,
    function(metric) {
      
      x <- sample_qc_star[[metric]]
      
      q1 <- quantile(x, 0.25, na.rm = TRUE)
      q3 <- quantile(x, 0.75, na.rm = TRUE)
      iqr <- IQR(x, na.rm = TRUE)
      
      data.frame(
        metric = metric,
        q1 = q1,
        q3 = q3,
        lower = q1 - 1.5 * iqr,
        upper = q3 + 1.5 * iqr,
        stringsAsFactors = FALSE
      )
    }
  )
)

star_iqr_thresholds
ga_lower <- star_iqr_thresholds$lower[
  star_iqr_thresholds$metric ==
    "gene_assignment_fraction"
]

unmapped_upper <- star_iqr_thresholds$upper[
  star_iqr_thresholds$metric ==
    "unmapped_fraction"
]

multimapping_upper <- star_iqr_thresholds$upper[
  star_iqr_thresholds$metric ==
    "multimapping_fraction"
]

noFeature_upper <- star_iqr_thresholds$upper[
  star_iqr_thresholds$metric ==
    "noFeature_fraction"
]

ambiguous_upper <- star_iqr_thresholds$upper[
  star_iqr_thresholds$metric ==
    "ambiguous_fraction"
]
sample_qc_star$low_assignment_flag <-
  sample_qc_star$gene_assignment_fraction <
  ga_lower

sample_qc_star$high_unmapped_flag <-
  sample_qc_star$unmapped_fraction >
  unmapped_upper

sample_qc_star$high_multimapping_flag <-
  sample_qc_star$multimapping_fraction >
  multimapping_upper

sample_qc_star$high_noFeature_flag <-
  sample_qc_star$noFeature_fraction >
  noFeature_upper

sample_qc_star$high_ambiguous_flag <-
  sample_qc_star$ambiguous_fraction >
  ambiguous_upper
sample_qc_star$n_star_flags <-
  rowSums(
    sample_qc_star[
      ,
      c(
        "low_assignment_flag",
        "high_unmapped_flag",
        "high_multimapping_flag",
        "high_noFeature_flag",
        "high_ambiguous_flag"
      )
    ]
  )
sample_qc_star$star_qc_flag <-
  sample_qc_star$n_star_flags > 0
table(sample_qc_star$low_assignment_flag)


table(sample_qc_star$high_unmapped_flag)

table(sample_qc_star$high_multimapping_flag)

table(sample_qc_star$high_noFeature_flag)

table(sample_qc_star$high_ambiguous_flag)

table(sample_qc_star$n_star_flags)

table(sample_qc_star$star_qc_flag)
table(
  raw_qc_flag = sample_qc_star$raw_qc_flag,
  star_qc_flag = sample_qc_star$star_qc_flag
)
convergent_qc_candidates <- sample_qc_star[
  sample_qc_star$raw_qc_flag &
    sample_qc_star$star_qc_flag,
  c(
    "patient_id",
    "assigned_gene_counts",
    "detected_genes",
    "zero_fraction",
    "gene_assignment_fraction",
    "unmapped_fraction",
    "multimapping_fraction",
    "noFeature_fraction",
    "ambiguous_fraction",
    "library_flag",
    "detected_flag",
    "low_assignment_flag",
    "high_unmapped_flag",
    "high_multimapping_flag",
    "high_noFeature_flag",
    "high_ambiguous_flag",
    "n_star_flags"
  )
]
convergent_qc_candidates <- convergent_qc_candidates[
  order(
    -convergent_qc_candidates$n_star_flags,
    convergent_qc_candidates$gene_assignment_fraction
  ),
]

convergent_qc_candidates
star_only_candidates <- sample_qc_star[
  !sample_qc_star$raw_qc_flag &
    sample_qc_star$star_qc_flag,
  c(
    "patient_id",
    "assigned_gene_counts",
    "detected_genes",
    "zero_fraction",
    "gene_assignment_fraction",
    "unmapped_fraction",
    "multimapping_fraction",
    "noFeature_fraction",
    "ambiguous_fraction",
    "n_star_flags"
  )
]

head(
  star_only_candidates[
    order(
      -star_only_candidates$n_star_flags,
      star_only_candidates$gene_assignment_fraction
    ),
  ],
  20
)
p_assignment <- ggplot(
  sample_qc_star,
  aes(x = gene_assignment_fraction)
) +
  geom_histogram(
    bins = 50,
    color = "white"
  ) +
  geom_vline(
    xintercept = ga_lower,
    linetype = "dashed"
  ) +
  labs(
    title = "Distribution of gene-assignment fraction",
    x = "Gene-assignment fraction",
    y = "Number of samples"
  ) +
  theme_minimal()

p_assignment
p_detected_assignment <- ggplot(
  sample_qc_star,
  aes(
    x = gene_assignment_fraction,
    y = detected_genes
  )
) +
  geom_point(
    alpha = 0.55
  ) +
  labs(
    title = "Detected genes versus gene-assignment fraction",
    x = "Gene-assignment fraction",
    y = "Detected genes"
  ) +
  theme_minimal()

p_detected_assignment
# ============================================================
# STEP 4E — Inspect available biospecimen metadata
# ============================================================

dim(biospecimen_raw)

names(biospecimen_raw)
grep(
  "ffpe|rin|rna|portion|analyte|plate|center|concentration|quantity|volume|preservation",
  names(biospecimen_raw),
  value = TRUE,
  ignore.case = TRUE
)
# ============================================================
# STEP 4E.1 — Inspect nested biospecimen structure
# ============================================================

class(biospecimen_raw$portions)

typeof(biospecimen_raw$portions)

str(
  biospecimen_raw$portions,
  max.level = 3
)
biospecimen_raw$portions[1:3]
sum(
  !is.na(biospecimen_raw$portions)
)
head(
  biospecimen_raw[
    ,
    c(
      "submitter_id",
      "sample_id",
      "sample_type",
      "specimen_type",
      "tissue_type",
      "preservation_method"
    )
  ],
  10
)
head(
  biospecimen_raw$submitter_id,
  20
)
biospecimen_raw$patient_id <-
  substr(
    biospecimen_raw$submitter_id,
    1,
    12
  )
head(
  biospecimen_raw[
    ,
    c(
      "patient_id",
      "submitter_id"
    )
  ]
)
table(
  biospecimen_raw$preservation_method,
  useNA = "ifany"
)
table(
  biospecimen_raw$sample_type,
  useNA = "ifany"
)
candidate_ids <-
  convergent_qc_candidates$patient_id

candidate_biospecimen <- biospecimen_raw[
  substr(
    biospecimen_raw$submitter_id,
    1,
    12
  ) %in% candidate_ids,
]

candidate_biospecimen[
  ,
  c(
    "submitter_id",
    "sample_type",
    "specimen_type",
    "tissue_type",
    "preservation_method",
    "freezing_method",
    "time_between_excision_and_freezing",
    "time_between_clamping_and_freezing"
  )
]
class(biospecimen_raw$portions)

typeof(biospecimen_raw$portions)

str(
  biospecimen_raw$portions,
  max.level = 3
)

head(
  biospecimen_raw$submitter_id,
  20
)
# ============================================================
# STEP 4E.2 — Map QC candidates to exact retained biospecimens
# ============================================================

candidate_manifest <- final_rna_manifest[
  final_rna_manifest$patient_id %in%
    convergent_qc_candidates$patient_id,
  c(
    "patient_id",
    "sample.submitter_id",
    "aliquot_id",
    "file_id",
    "file_name"
  )
]

candidate_manifest <- candidate_manifest[
  match(
    convergent_qc_candidates$patient_id,
    candidate_manifest$patient_id
  ),
]

nrow(candidate_manifest)

anyDuplicated(candidate_manifest$patient_id)

anyNA(candidate_manifest$sample.submitter_id)

anyNA(candidate_manifest$aliquot_id)
table(
  biospecimen_raw$preservation_method,
  useNA = "ifany"
)
candidate_sample_match <- match(
  candidate_manifest$sample.submitter_id,
  biospecimen_raw$submitter_id
)

sum(is.na(candidate_sample_match))
candidate_sample_metadata <- data.frame(
  patient_id =
    candidate_manifest$patient_id,
  
  retained_sample_id =
    candidate_manifest$sample.submitter_id,
  
  retained_aliquot_id =
    candidate_manifest$aliquot_id,
  
  sample_type =
    biospecimen_raw$sample_type[
      candidate_sample_match
    ],
  
  specimen_type =
    biospecimen_raw$specimen_type[
      candidate_sample_match
    ],
  
  tissue_type =
    biospecimen_raw$tissue_type[
      candidate_sample_match
    ],
  
  preservation_method =
    biospecimen_raw$preservation_method[
      candidate_sample_match
    ],
  
  freezing_method =
    biospecimen_raw$freezing_method[
      candidate_sample_match
    ],
  
  time_excision_to_freezing =
    biospecimen_raw$time_between_excision_and_freezing[
      candidate_sample_match
    ],
  
  time_clamping_to_freezing =
    biospecimen_raw$time_between_clamping_and_freezing[
      candidate_sample_match
    ],
  
  stringsAsFactors = FALSE
)
candidate_sample_metadata
table(
  candidate_sample_metadata$sample_type,
  useNA = "ifany"
)

table(
  candidate_sample_metadata$specimen_type,
  useNA = "ifany"
)

table(
  candidate_sample_metadata$tissue_type,
  useNA = "ifany"
)

table(
  candidate_sample_metadata$preservation_method,
  useNA = "ifany"
)
a3od_row <- match(
  candidate_manifest$sample.submitter_id[
    candidate_manifest$patient_id ==
      "TCGA-AC-A3OD"
  ],
  biospecimen_raw$submitter_id
)

a3od_row
a3od_portions <-
  biospecimen_raw$portions[[a3od_row]]

a3od_portions[
  ,
  setdiff(
    names(a3od_portions),
    "analytes"
  ),
  drop = FALSE
]
nrow(a3od_portions)

names(a3od_portions)
str(
  a3od_portions$analytes,
  max.level = 4
)
str(
  a3od_portions$analytes,
  max.level = 4
)
lapply(
  a3od_portions$analytes,
  names
)
which(
  vapply(
    a3od_portions$analytes,
    is.data.frame,
    logical(1)
  )
)
a3od_analyte_index <- which(
  vapply(
    a3od_portions$analytes,
    is.data.frame,
    logical(1)
  )
)[1]

a3od_analytes <-
  a3od_portions$analytes[[a3od_analyte_index]]

names(a3od_analytes)

a3od_analytes
nrow(candidate_manifest)
sum(is.na(candidate_sample_match))

table(
  candidate_sample_metadata$preservation_method,
  useNA = "ifany"
)

table(
  candidate_sample_metadata$specimen_type,
  useNA = "ifany"
)

a3od_portions[
  ,
  setdiff(
    names(a3od_portions),
    "analytes"
  ),
  drop = FALSE
]

names(a3od_analytes)

a3od_analytes

grep(
  "aliquot",
  names(a3od_analytes),
  value = TRUE,
  ignore.case = TRUE
)
names(a3od_analytes)
a3od_rna <- a3od_analytes[
  a3od_analytes$analyte_type == "RNA",
  ,
  drop = FALSE
]

a3od_rna[
  ,
  setdiff(
    names(a3od_rna),
    "aliquots"
  ),
  drop = FALSE
]
nrow(a3od_rna)
a3od_rna_aliquots <-
  a3od_rna$aliquots[[1]]

a3od_rna_aliquots
a3od_manifest <- final_rna_manifest[
  final_rna_manifest$patient_id == "TCGA-AC-A3OD",
  c(
    "patient_id",
    "sample.submitter_id",
    "aliquot_id",
    "file_id",
    "file_name"
  )
]

a3od_manifest
a3od_aliquot_match <- match(
  a3od_manifest$aliquot_id,
  a3od_rna_aliquots$submitter_id
)

a3od_aliquot_match
a3od_retained_aliquot <-
  a3od_rna_aliquots[
    a3od_aliquot_match,
    ,
    drop = FALSE
  ]

a3od_retained_aliquot
# ============================================================
# STEP 4E.3 — Extract exact biospecimen metadata
# for convergent QC candidates
# ============================================================

extract_rna_biospecimen <- function(patient_id,
                                    sample_id,
                                    aliquot_barcode,
                                    biospecimen_raw) {
  
  # ----------------------------------------------------------
  # 1. Find exact sample/vial
  # ----------------------------------------------------------
  
  sample_index <- match(
    sample_id,
    biospecimen_raw$submitter_id
  )
  
  if (is.na(sample_index)) {
    return(
      data.frame(
        patient_id = patient_id,
        sample_id = sample_id,
        aliquot_id = aliquot_barcode,
        mapping_status = "sample_not_found",
        stringsAsFactors = FALSE
      )
    )
  }
  
  portions <- biospecimen_raw$portions[[sample_index]]
  
  if (
    is.null(portions) ||
    !is.data.frame(portions) ||
    !"analytes" %in% names(portions)
  ) {
    return(
      data.frame(
        patient_id = patient_id,
        sample_id = sample_id,
        aliquot_id = aliquot_barcode,
        mapping_status = "no_nested_portion_data",
        stringsAsFactors = FALSE
      )
    )
  }
  
  # ----------------------------------------------------------
  # 2. Search all portions → analytes → aliquots
  #    for the exact retained RNA aliquot
  # ----------------------------------------------------------
  
  for (p in seq_len(nrow(portions))) {
    
    analytes <- portions$analytes[[p]]
    
    if (
      is.null(analytes) ||
      !is.data.frame(analytes) ||
      !"aliquots" %in% names(analytes)
    ) {
      next
    }
    
    for (a in seq_len(nrow(analytes))) {
      
      aliquots <- analytes$aliquots[[a]]
      
      if (
        is.null(aliquots) ||
        !is.data.frame(aliquots) ||
        !"submitter_id" %in% names(aliquots)
      ) {
        next
      }
      
      q <- match(
        aliquot_barcode,
        aliquots$submitter_id
      )
      
      if (!is.na(q)) {
        
        # Safe extractor for columns that may not exist
        get_value <- function(df, column, row = 1) {
          if (column %in% names(df)) {
            df[[column]][row]
          } else {
            NA
          }
        }
        
        return(
          data.frame(
            patient_id = patient_id,
            
            sample_id = sample_id,
            
            sample_type =
              biospecimen_raw$sample_type[sample_index],
            
            specimen_type =
              biospecimen_raw$specimen_type[sample_index],
            
            preservation_method =
              biospecimen_raw$preservation_method[sample_index],
            
            portion_id =
              get_value(portions, "submitter_id", p),
            
            portion_is_ffpe =
              get_value(portions, "is_ffpe", p),
            
            portion_weight =
              get_value(portions, "weight", p),
            
            analyte_id =
              get_value(analytes, "submitter_id", a),
            
            analyte_type =
              get_value(analytes, "analyte_type", a),
            
            experimental_protocol =
              get_value(
                analytes,
                "experimental_protocol_type",
                a
              ),
            
            rna_integrity_number =
              get_value(
                analytes,
                "rna_integrity_number",
                a
              ),
            
            ribosomal_28s_18s_ratio =
              get_value(
                analytes,
                "ribosomal_rna_28s_18s_ratio",
                a
              ),
            
            a260_a280_ratio =
              get_value(
                analytes,
                "a260_a280_ratio",
                a
              ),
            
            rna_concentration =
              get_value(
                analytes,
                "concentration",
                a
              ),
            
            aliquot_id =
              aliquots$submitter_id[q],
            
            aliquot_source_center =
              get_value(
                aliquots,
                "source_center",
                q
              ),
            
            aliquot_concentration =
              get_value(
                aliquots,
                "concentration",
                q
              ),
            
            aliquot_quantity =
              get_value(
                aliquots,
                "aliquot_quantity",
                q
              ),
            
            aliquot_volume =
              get_value(
                aliquots,
                "aliquot_volume",
                q
              ),
            
            mapping_status = "exact_aliquot_match",
            
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  
  # exact aliquot was not found
  data.frame(
    patient_id = patient_id,
    sample_id = sample_id,
    aliquot_id = aliquot_barcode,
    mapping_status = "aliquot_not_found",
    stringsAsFactors = FALSE
  )
}
candidate_biospecimen_qc <- do.call(
  rbind,
  lapply(
    seq_len(nrow(candidate_manifest)),
    function(i) {
      
      extract_rna_biospecimen(
        patient_id =
          candidate_manifest$patient_id[i],
        
        sample_id =
          candidate_manifest$sample.submitter_id[i],
        
        aliquot_barcode =
          candidate_manifest$aliquot_id[i],
        
        biospecimen_raw =
          biospecimen_raw
      )
    }
  )
)
dim(candidate_biospecimen_qc)

table(
  candidate_biospecimen_qc$mapping_status,
  useNA = "ifany"
)
candidate_qc_integrated <- merge(
  convergent_qc_candidates,
  candidate_biospecimen_qc,
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)
candidate_qc_integrated <-
  candidate_qc_integrated[
    match(
      convergent_qc_candidates$patient_id,
      candidate_qc_integrated$patient_id
    ),
  ]
nrow(candidate_qc_integrated)

identical(
  candidate_qc_integrated$patient_id,
  convergent_qc_candidates$patient_id
)
table(
  candidate_qc_integrated$portion_is_ffpe,
  useNA = "ifany"
)
summary(
  candidate_qc_integrated$rna_integrity_number
)
candidate_qc_integrated[
  order(
    candidate_qc_integrated$rna_integrity_number,
    na.last = TRUE
  ),
  c(
    "patient_id",
    "n_star_flags",
    "detected_genes",
    "gene_assignment_fraction",
    "portion_is_ffpe",
    "rna_integrity_number",
    "ribosomal_28s_18s_ratio",
    "a260_a280_ratio",
    "aliquot_source_center"
  )
]
severe_ids <- c(
  "TCGA-AC-A3OD",
  "TCGA-A7-A13G",
  "TCGA-A7-A26F",
  "TCGA-A7-A26I",
  "TCGA-AC-A3QQ",
  "TCGA-AC-A2QH"
)

severe_biospecimen_qc <-
  candidate_qc_integrated[
    candidate_qc_integrated$patient_id %in%
      severe_ids,
    c(
      "patient_id",
      "n_star_flags",
      "assigned_gene_counts",
      "detected_genes",
      "gene_assignment_fraction",
      "unmapped_fraction",
      "multimapping_fraction",
      "noFeature_fraction",
      "ambiguous_fraction",
      "preservation_method",
      "portion_id",
      "portion_is_ffpe",
      "analyte_id",
      "rna_integrity_number",
      "ribosomal_28s_18s_ratio",
      "a260_a280_ratio",
      "aliquot_id",
      "aliquot_source_center"
    )
  ]

severe_biospecimen_qc
table(
  candidate_qc_integrated$aliquot_source_center,
  useNA = "ifany"
)
table(
  severe_biospecimen_qc$aliquot_source_center,
  useNA = "ifany"
)
table(
  candidate_biospecimen_qc$mapping_status,
  useNA = "ifany"
)

table(
  candidate_qc_integrated$portion_is_ffpe,
  useNA = "ifany"
)

summary(
  candidate_qc_integrated$rna_integrity_number
)

severe_biospecimen_qc
# ============================================================
# STEP 5 — EXPLORATORY PRE-QC FILTERING (not the final inference input)
# STEP 5.1 — Create edgeR analysis object
# ============================================================
if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop("edgeR is required; install it separately before running this section.")
}
library(edgeR)
packageVersion("edgeR")
dge <- DGEList(
  counts = raw_count_matrix
)

dim(dge)
identical(
  colnames(dge$counts),
  matched_cohort$patient_id
)

all(
  dge$counts == raw_count_matrix
)
all_zero_genes <- rowSums(
  raw_count_matrix
) == 0

sum(all_zero_genes)
genes_detected_samples <- rowSums(
  raw_count_matrix > 0
)

summary(
  genes_detected_samples
)
raw_cpm <- cpm(
  dge,
  normalized.lib.sizes = FALSE
)

dim(raw_cpm)

anyNA(raw_cpm)
samples_cpm_ge_1 <- rowSums(
  raw_cpm >= 1
)

summary(
  samples_cpm_ge_1
)
filter_sensitivity <- data.frame(
  
  criterion = c(
    "CPM >= 1 in >= 5 samples",
    "CPM >= 1 in >= 10 samples",
    "CPM >= 1 in >= 20 samples",
    "CPM >= 1 in >= 50 samples",
    "CPM >= 1 in >= 100 samples"
  ),
  
  genes_retained = c(
    sum(samples_cpm_ge_1 >= 5),
    sum(samples_cpm_ge_1 >= 10),
    sum(samples_cpm_ge_1 >= 20),
    sum(samples_cpm_ge_1 >= 50),
    sum(samples_cpm_ge_1 >= 100)
  )
)

filter_sensitivity$genes_removed <-
  nrow(raw_count_matrix) -
  filter_sensitivity$genes_retained

filter_sensitivity$retained_percent <-
  100 *
  filter_sensitivity$genes_retained /
  nrow(raw_count_matrix)

filter_sensitivity
# ============================================================
# STEP 5.4 — Design-aware low-expression filtering
# ============================================================

stage_idx <- !is.na(
  matched_cohort$stage_group
)

sum(stage_idx)
stage_meta <- droplevels(
  matched_cohort[stage_idx, ]
)

table(stage_meta$stage_group)
design_stage <- model.matrix(
  ~ stage_group,
  data = stage_meta
)

dim(design_stage)
dge_stage <- dge[
  ,
  stage_idx,
  keep.lib.sizes = FALSE
]

dim(dge_stage)
keep_filterByExpr_stage <- filterByExpr(
  dge_stage,
  design = design_stage
)

table(keep_filterByExpr_stage)

sum(keep_filterByExpr_stage)
manual_keep_20 <-
  samples_cpm_ge_1 >= 20
table(
  CPM_1_in_20 = manual_keep_20,
  filterByExpr = keep_filterByExpr_stage
)
c(
  manual_CPM20 =
    sum(manual_keep_20),
  
  filterByExpr_stage =
    sum(keep_filterByExpr_stage)
)
sum(
  manual_keep_20 &
    keep_filterByExpr_stage
)
sum(
  manual_keep_20 &
    !keep_filterByExpr_stage
)

sum(
  !manual_keep_20 &
    keep_filterByExpr_stage
)
age_idx <- !is.na(
  matched_cohort$age_years
)

sum(age_idx)
age_meta <- matched_cohort[
  age_idx,
]

design_age <- model.matrix(
  ~ age_years,
  data = age_meta
)

dim(design_age)
dge_age <- dge[
  ,
  age_idx,
  keep.lib.sizes = FALSE
]
keep_filterByExpr_age <- filterByExpr(
  dge_age,
  design = design_age
)

table(keep_filterByExpr_age)

sum(keep_filterByExpr_age)
table(
  Age = keep_filterByExpr_age,
  Stage = keep_filterByExpr_stage
)
c(
  age = sum(keep_filterByExpr_age),
  stage = sum(keep_filterByExpr_stage),
  both = sum(
    keep_filterByExpr_age &
      keep_filterByExpr_stage
  ),
  either = sum(
    keep_filterByExpr_age |
      keep_filterByExpr_stage
  )
)
# ============================================================
# STEP 5.6 — Create filtered analysis-specific objects
# ============================================================

dge_age_filtered <- dge_age[
  keep_filterByExpr_age,
  ,
  keep.lib.sizes = FALSE
]

dim(dge_age_filtered)
dge_stage_filtered <- dge_stage[
  keep_filterByExpr_stage,
  ,
  keep.lib.sizes = FALSE
]

dim(dge_stage_filtered)
identical(
  colnames(dge_age_filtered$counts),
  age_meta$patient_id
)

identical(
  colnames(dge_stage_filtered$counts),
  stage_meta$patient_id
)
c(
  age_genes = nrow(dge_age_filtered),
  stage_genes = nrow(dge_stage_filtered)
)
# ============================================================
# STEP 6.1 — TMM normalization factors
# ============================================================

dge_age_norm <- calcNormFactors(
  dge_age_filtered,
  method = "TMM"
)

dge_age_norm$samples
dge_stage_norm <- calcNormFactors(
  dge_stage_filtered,
  method = "TMM"
)

dge_stage_norm$samples
summary(
  dge_age_norm$samples$norm.factors
)

summary(
  dge_stage_norm$samples$norm.factors
)
range(
  dge_age_norm$samples$norm.factors
)

range(
  dge_stage_norm$samples$norm.factors
)
prod(
  dge_age_norm$samples$norm.factors
)

prod(
  dge_stage_norm$samples$norm.factors
)
# ============================================================
# STEP 6.2 — Effective library sizes after TMM
# ============================================================

dge_age_norm$samples$effective_lib_size <-
  dge_age_norm$samples$lib.size *
  dge_age_norm$samples$norm.factors

dge_stage_norm$samples$effective_lib_size <-
  dge_stage_norm$samples$lib.size *
  dge_stage_norm$samples$norm.factors
summary(
  dge_age_norm$samples$effective_lib_size
)

summary(
  dge_stage_norm$samples$effective_lib_size
)
cor(
  dge_age_norm$samples$lib.size,
  dge_age_norm$samples$effective_lib_size,
  method = "spearman"
)

cor(
  dge_stage_norm$samples$lib.size,
  dge_stage_norm$samples$effective_lib_size,
  method = "spearman"
)
head(
  dge_age_norm$samples[
    order(dge_age_norm$samples$norm.factors),
  ],
  10
)

tail(
  dge_age_norm$samples[
    order(dge_age_norm$samples$norm.factors),
  ],
  10
)
head(
  dge_stage_norm$samples[
    order(dge_stage_norm$samples$norm.factors),
  ],
  10
)

tail(
  dge_stage_norm$samples[
    order(dge_stage_norm$samples$norm.factors),
  ],
  10
)
library(limma)
# ============================================================
# STEP 6.3 — voom transformation
# ============================================================

voom_age <- voom(
  dge_age_norm,
  design = design_age,
  plot = TRUE
)
voom_stage <- voom(
  dge_stage_norm,
  design = design_stage,
  plot = TRUE
)
class(voom_age)
dim(voom_age$E)
dim(voom_age$weights)

class(voom_stage)
dim(voom_stage$E)
dim(voom_stage$weights)
identical(
  colnames(voom_age$E),
  age_meta$patient_id
)

identical(
  colnames(voom_stage$E),
  stage_meta$patient_id
)
summary(as.vector(voom_age$weights))

summary(as.vector(voom_stage$weights))
# ============================================================
# STEP 7.1 — PCA on voom-transformed expression
# ============================================================

pipeline_event(ctx, "diagnostic_pca_age")
pca_age <- prcomp(
  t(voom_age$E),
  center = TRUE,
  scale. = FALSE
)

pipeline_event(ctx, "diagnostic_pca_stage")
pca_stage <- prcomp(
  t(voom_stage$E),
  center = TRUE,
  scale. = FALSE
)
age_var_explained <-
  100 * pca_age$sdev^2 /
  sum(pca_age$sdev^2)

stage_var_explained <-
  100 * pca_stage$sdev^2 /
  sum(pca_stage$sdev^2)

age_var_explained[1:10]

stage_var_explained[1:10]
pca_age_df <- data.frame(
  patient_id = rownames(pca_age$x),
  PC1 = pca_age$x[, 1],
  PC2 = pca_age$x[, 2],
  PC3 = pca_age$x[, 3],
  stringsAsFactors = FALSE
)

pca_stage_df <- data.frame(
  patient_id = rownames(pca_stage$x),
  PC1 = pca_stage$x[, 1],
  PC2 = pca_stage$x[, 2],
  PC3 = pca_stage$x[, 3],
  stringsAsFactors = FALSE
)
pca_age_df <- merge(
  pca_age_df,
  age_meta,
  by = "patient_id",
  sort = FALSE
)

pca_stage_df <- merge(
  pca_stage_df,
  stage_meta,
  by = "patient_id",
  sort = FALSE
)
pca_age_df <- pca_age_df[
  match(
    rownames(pca_age$x),
    pca_age_df$patient_id
  ),
]

pca_stage_df <- pca_stage_df[
  match(
    rownames(pca_stage$x),
    pca_stage_df$patient_id
  ),
]
identical(
  pca_age_df$patient_id,
  rownames(pca_age$x)
)

identical(
  pca_stage_df$patient_id,
  rownames(pca_stage$x)
)
qc_candidates <- c(
  "TCGA-AC-A3QQ",
  "TCGA-AC-A3OD",
  "TCGA-A7-A13G",
  "TCGA-A7-A26F",
  "TCGA-A7-A26I",
  "TCGA-AC-A2QH"
)

pca_age_df$qc_candidate <-
  pca_age_df$patient_id %in% qc_candidates

pca_stage_df$qc_candidate <-
  pca_stage_df$patient_id %in% qc_candidates
table(pca_age_df$qc_candidate)

table(pca_stage_df$qc_candidate)
dim(pca_age$x)

dim(pca_stage$x)
# ============================================================
# STEP 7.2 — PCA sample-level QC visualization
# ============================================================

library(ggplot2)

p_age_qc <- ggplot(
  pca_age_df,
  aes(x = PC1, y = PC2)
) +
  geom_point(
    alpha = 0.55,
    size = 1.7
  ) +
  labs(
    title = "TCGA-BRCA — Age Analysis",
    subtitle = "Post-normalization PCA for sample-level QC",
    x = paste0(
      "PC1 (",
      round(age_var_explained[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(age_var_explained[2], 1),
      "%)"
    )
  ) +
  theme_classic()

p_age_qc
p_stage_qc <- ggplot(
  pca_stage_df,
  aes(x = PC1, y = PC2)
) +
  geom_point(
    alpha = 0.55,
    size = 1.7
  ) +
  labs(
    title = "TCGA-BRCA — Stage Analysis",
    subtitle = "Post-normalization PCA for sample-level QC",
    x = paste0(
      "PC1 (",
      round(stage_var_explained[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(stage_var_explained[2], 1),
      "%)"
    )
  ) +
  theme_classic()

p_stage_qc
p_age_candidates <- ggplot(
  pca_age_df,
  aes(x = PC1, y = PC2)
) +
  geom_point(
    alpha = 0.35,
    size = 1.6
  ) +
  geom_point(
    data = subset(
      pca_age_df,
      qc_candidate
    ),
    size = 3
  ) +
  geom_text(
    data = subset(
      pca_age_df,
      qc_candidate
    ),
    aes(label = patient_id),
    hjust = -0.1,
    vjust = -0.5,
    size = 3
  ) +
  labs(
    title = "TCGA-BRCA — Age Analysis",
    subtitle = "Previously identified QC candidates",
    x = paste0(
      "PC1 (",
      round(age_var_explained[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(age_var_explained[2], 1),
      "%)"
    )
  ) +
  theme_classic()

p_age_candidates
p_stage_candidates <- ggplot(
  pca_stage_df,
  aes(x = PC1, y = PC2)
) +
  geom_point(
    alpha = 0.35,
    size = 1.6
  ) +
  geom_point(
    data = subset(
      pca_stage_df,
      qc_candidate
    ),
    size = 3
  ) +
  geom_text(
    data = subset(
      pca_stage_df,
      qc_candidate
    ),
    aes(label = patient_id),
    hjust = -0.1,
    vjust = -0.5,
    size = 3
  ) +
  labs(
    title = "TCGA-BRCA — Stage Analysis",
    subtitle = "Previously identified QC candidates",
    x = paste0(
      "PC1 (",
      round(stage_var_explained[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(stage_var_explained[2], 1),
      "%)"
    )
  ) +
  theme_classic()

p_stage_candidates
# ============================================================
# STEP 7.3 — Quantitative multivariate PCA screening
# ============================================================

age_cumulative_variance <- cumsum(age_var_explained)
stage_cumulative_variance <- cumsum(stage_var_explained)

age_npc_50 <- which(
  age_cumulative_variance >= 50
)[1]

stage_npc_50 <- which(
  stage_cumulative_variance >= 50
)[1]

age_npc_50
stage_npc_50

age_cumulative_variance[age_npc_50]
stage_cumulative_variance[stage_npc_50]
age_scores <- pca_age$x[
  ,
  1:age_npc_50,
  drop = FALSE
]

stage_scores <- pca_stage$x[
  ,
  1:stage_npc_50,
  drop = FALSE
]


age_scores_z <- scale(age_scores)

stage_scores_z <- scale(stage_scores)


age_pca_distance <- sqrt(
  rowSums(age_scores_z^2)
)

stage_pca_distance <- sqrt(
  rowSums(stage_scores_z^2)
)
pca_age_df$pca_distance <-
  age_pca_distance[
    match(
      pca_age_df$patient_id,
      names(age_pca_distance)
    )
  ]

pca_stage_df$pca_distance <-
  stage_pca_distance[
    match(
      pca_stage_df$patient_id,
      names(stage_pca_distance)
    )
  ]
anyNA(pca_age_df$pca_distance)
anyNA(pca_stage_df$pca_distance)

summary(pca_age_df$pca_distance)
summary(pca_stage_df$pca_distance)
age_candidate_rank <- pca_age_df[
  pca_age_df$qc_candidate,
  c(
    "patient_id",
    "PC1",
    "PC2",
    "PC3",
    "pca_distance"
  )
]

age_candidate_rank$rank <-
  rank(
    -pca_age_df$pca_distance,
    ties.method = "min"
  )[
    pca_age_df$qc_candidate
  ]

age_candidate_rank <- age_candidate_rank[
  order(age_candidate_rank$rank),
]

age_candidate_rank
stage_candidate_rank <- pca_stage_df[
  pca_stage_df$qc_candidate,
  c(
    "patient_id",
    "PC1",
    "PC2",
    "PC3",
    "pca_distance"
  )
]

stage_candidate_rank$rank <-
  rank(
    -pca_stage_df$pca_distance,
    ties.method = "min"
  )[
    pca_stage_df$qc_candidate
  ]

stage_candidate_rank <- stage_candidate_rank[
  order(stage_candidate_rank$rank),
]

stage_candidate_rank
head(
  pca_age_df[
    order(
      pca_age_df$pca_distance,
      decreasing = TRUE
    ),
    c(
      "patient_id",
      "pca_distance"
    )
  ],
  20
)
head(
  pca_stage_df[
    order(
      pca_stage_df$pca_distance,
      decreasing = TRUE
    ),
    c(
      "patient_id",
      "pca_distance"
    )
  ],
  20
)
# ============================================================
pipeline_event(ctx, "canonical_preprocessing")
# STEP 8 — CANONICAL final Age/Stage cohorts and preprocessing after QC
# ============================================================

final_qc_exclusions <- c(
  "TCGA-AC-A2QH",
  "TCGA-AC-A3OD",
  "TCGA-AC-A3QQ",
  "TCGA-A7-A13G",
  "TCGA-A7-A26F",
  "TCGA-A7-A26I"
)
# One post-QC parent cohort, in explicit metadata patient order.
stopifnot(
  is.matrix(raw_count_matrix), is.numeric(raw_count_matrix),
  nrow(raw_count_matrix) > 0L, ncol(raw_count_matrix) > 0L,
  !anyNA(raw_count_matrix), all(is.finite(raw_count_matrix)),
  all(raw_count_matrix >= 0), all(raw_count_matrix == round(raw_count_matrix)),
  !is.null(rownames(raw_count_matrix)), !is.null(colnames(raw_count_matrix)),
  !anyNA(rownames(raw_count_matrix)), !anyNA(colnames(raw_count_matrix)),
  all(nzchar(rownames(raw_count_matrix))), all(nzchar(colnames(raw_count_matrix))),
  !anyDuplicated(rownames(raw_count_matrix)), !anyDuplicated(colnames(raw_count_matrix)),
  is.character(matched_cohort$patient_id), !anyNA(matched_cohort$patient_id),
  all(nzchar(matched_cohort$patient_id)), !anyDuplicated(matched_cohort$patient_id),
  setequal(colnames(raw_count_matrix), matched_cohort$patient_id),
  all(final_qc_exclusions %in% matched_cohort$patient_id)
)
matched_cohort_final <- matched_cohort[
  !matched_cohort$patient_id %in% final_qc_exclusions, , drop = FALSE
]
raw_count_matrix_final <- raw_count_matrix[
  , matched_cohort_final$patient_id, drop = FALSE
]
stopifnot(
  identical(colnames(raw_count_matrix_final), matched_cohort_final$patient_id),
  !any(final_qc_exclusions %in% colnames(raw_count_matrix_final))
)

# Age eligibility is age availability; Stage availability is not required.
age_meta_final <- matched_cohort_final[
  !is.na(matched_cohort_final$age_years), , drop = FALSE
]
# Stage eligibility is harmonized Stage availability; age availability is not required.
stage_meta_final <- matched_cohort_final[
  !is.na(matched_cohort_final$stage_group), , drop = FALSE
]
stage_meta_final$stage_group <- factor(
  as.character(stage_meta_final$stage_group),
  levels = c("Stage I", "Stage II", "Stage III", "Stage IV")
)
rownames(age_meta_final) <- age_meta_final$patient_id
rownames(stage_meta_final) <- stage_meta_final$patient_id
stopifnot(
  nrow(age_meta_final) > 2L, nrow(stage_meta_final) > 4L,
  all(is.finite(age_meta_final$age_years)), !anyNA(stage_meta_final$stage_group)
)

# Select by patient_id, never by incidental row names or logical column position.
raw_count_age_final <- raw_count_matrix_final[
  , age_meta_final$patient_id, drop = FALSE
]
raw_count_stage_final <- raw_count_matrix_final[
  , stage_meta_final$patient_id, drop = FALSE
]
dge_age_final <- DGEList(counts = raw_count_age_final)
dge_stage_final <- DGEList(counts = raw_count_stage_final)
design_age_final <- model.matrix(~ age_years, data = age_meta_final)
design_stage_final <- model.matrix(
  ~ stage_group, data = stage_meta_final,
  contrasts.arg = list(stage_group = contr.treatment(
    c("Stage I", "Stage II", "Stage III", "Stage IV"), base = 1L
  ))
)
stopifnot(
  identical(rownames(design_age_final), colnames(raw_count_age_final)),
  identical(rownames(design_stage_final), colnames(raw_count_stage_final)),
  nrow(design_age_final) == nrow(age_meta_final),
  nrow(design_stage_final) == nrow(stage_meta_final),
  qr(design_age_final)$rank == ncol(design_age_final),
  qr(design_stage_final)$rank == ncol(design_stage_final)
)

# Filter, TMM and voom remain specific to each analysis cohort and design.
keep_age_final <- filterByExpr(dge_age_final, design = design_age_final)
keep_stage_final <- filterByExpr(dge_stage_final, design = design_stage_final)
stopifnot(any(keep_age_final), any(keep_stage_final))
dge_age_final_filtered <- dge_age_final[
  keep_age_final, , keep.lib.sizes = FALSE
]
dge_stage_final_filtered <- dge_stage_final[
  keep_stage_final, , keep.lib.sizes = FALSE
]
# ============================================================
# STEP 8.4 — Final TMM normalization
# ============================================================

dge_age_final_norm <- calcNormFactors(
  dge_age_final_filtered,
  method = "TMM"
)

dge_stage_final_norm <- calcNormFactors(
  dge_stage_final_filtered,
  method = "TMM"
)
summary(
  dge_age_final_norm$samples$norm.factors
)

summary(
  dge_stage_final_norm$samples$norm.factors
)

range(
  dge_age_final_norm$samples$norm.factors
)

range(
  dge_stage_final_norm$samples$norm.factors
)

prod(
  dge_age_final_norm$samples$norm.factors
)

prod(
  dge_stage_final_norm$samples$norm.factors
)
dge_age_final_norm$samples$effective_lib_size <-
  dge_age_final_norm$samples$lib.size *
  dge_age_final_norm$samples$norm.factors

dge_stage_final_norm$samples$effective_lib_size <-
  dge_stage_final_norm$samples$lib.size *
  dge_stage_final_norm$samples$norm.factors

summary(
  dge_age_final_norm$samples$effective_lib_size
)

summary(
  dge_stage_final_norm$samples$effective_lib_size
)
# ============================================================
# STEP 8.5 — Final voom transformation
# ============================================================

library(limma)

voom_age_final <- voom(
  dge_age_final_norm,
  design = design_age_final,
  plot = FALSE
)
voom_stage_final <- voom(
  dge_stage_final_norm,
  design = design_stage_final,
  plot = FALSE
)
class(voom_age_final)
class(voom_stage_final)

dim(voom_age_final$E)
dim(voom_age_final$weights)

dim(voom_stage_final$E)
dim(voom_stage_final$weights)
identical(
  colnames(voom_age_final$E),
  age_meta_final$patient_id
)

identical(
  colnames(voom_stage_final$E),
  stage_meta_final$patient_id
)

identical(
  colnames(voom_age_final$weights),
  age_meta_final$patient_id
)

identical(
  colnames(voom_stage_final$weights),
  stage_meta_final$patient_id
)
summary(
  as.vector(voom_age_final$weights)
)

summary(
  as.vector(voom_stage_final$weights)
)

anyNA(voom_age_final$E)
anyNA(voom_stage_final$E)

anyNA(voom_age_final$weights)
anyNA(voom_stage_final$weights)
ncol(voom_age_final$E) ==
  nrow(design_age_final)

ncol(voom_stage_final$E) ==
  nrow(design_stage_final)
identical(
  dim(voom_age_final$weights),
  dim(voom_age_final$E)
)

identical(
  dim(voom_stage_final$weights),
  dim(voom_stage_final$E)
)
stopifnot(
  identical(colnames(voom_age_final$E), age_meta_final$patient_id),
  identical(colnames(voom_stage_final$E), stage_meta_final$patient_id),
  identical(dim(voom_age_final$E), dim(voom_age_final$weights)),
  identical(dim(voom_stage_final$E), dim(voom_stage_final$weights)),
  identical(rownames(voom_age_final$E), rownames(dge_age_final_norm$counts)),
  identical(rownames(voom_stage_final$E), rownames(dge_stage_final_norm$counts)),
  all(is.finite(voom_age_final$E)), all(is.finite(voom_stage_final$E)),
  all(is.finite(voom_age_final$weights)), all(voom_age_final$weights > 0),
  all(is.finite(voom_stage_final$weights)), all(voom_stage_final$weights > 0),
  all(is.finite(dge_age_final_norm$samples$norm.factors)),
  all(dge_age_final_norm$samples$norm.factors > 0),
  all(is.finite(dge_stage_final_norm$samples$norm.factors)),
  all(dge_stage_final_norm$samples$norm.factors > 0)
)
# ============================================================
# STEP 8.6 — Final PCA sanity check
# ============================================================

pipeline_event(ctx, "diagnostic_pca_age_final")
pca_age_final <- prcomp(
  t(voom_age_final$E),
  center = TRUE,
  scale. = FALSE
)

pipeline_event(ctx, "diagnostic_pca_stage_final")
pca_stage_final <- prcomp(
  t(voom_stage_final$E),
  center = TRUE,
  scale. = FALSE
)
age_final_var_explained <-
  100 * pca_age_final$sdev^2 /
  sum(pca_age_final$sdev^2)

stage_final_var_explained <-
  100 * pca_stage_final$sdev^2 /
  sum(pca_stage_final$sdev^2)

age_final_var_explained[1:10]

stage_final_var_explained[1:10]

sum(age_final_var_explained[1:2])

sum(stage_final_var_explained[1:2])
pca_age_final_df <- data.frame(
  patient_id = rownames(pca_age_final$x),
  PC1 = pca_age_final$x[, 1],
  PC2 = pca_age_final$x[, 2],
  age_years = age_meta_final$age_years
)

pca_stage_final_df <- data.frame(
  patient_id = rownames(pca_stage_final$x),
  PC1 = pca_stage_final$x[, 1],
  PC2 = pca_stage_final$x[, 2],
  stage_group = stage_meta_final$stage_group
)
identical(
  pca_age_final_df$patient_id,
  age_meta_final$patient_id
)

identical(
  pca_stage_final_df$patient_id,
  stage_meta_final$patient_id
)
plot(
  pca_age_final_df$PC1,
  pca_age_final_df$PC2,
  pch = 16,
  cex = 0.6,
  xlab = paste0(
    "PC1 (",
    round(age_final_var_explained[1], 1),
    "%)"
  ),
  ylab = paste0(
    "PC2 (",
    round(age_final_var_explained[2], 1),
    "%)"
  ),
  main = "Final PCA — Age Analysis"
)
plot(
  pca_stage_final_df$PC1,
  pca_stage_final_df$PC2,
  pch = 16,
  cex = 0.6,
  xlab = paste0(
    "PC1 (",
    round(stage_final_var_explained[1], 1),
    "%)"
  ),
  ylab = paste0(
    "PC2 (",
    round(stage_final_var_explained[2], 1),
    "%)"
  ),
  main = "Final PCA — Stage Analysis"
)
# ============================================================
# Save final analysis-ready objects
# ============================================================

dir.create(
  file.path(output_dir, "data_processed/final_analysis"),
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  voom_age_final,
  file.path(output_dir, "data_processed/final_analysis/TCGA_BRCA_age_voom_final.rds")
)

saveRDS(
  voom_stage_final,
  file.path(output_dir, "data_processed/final_analysis/TCGA_BRCA_stage_voom_final.rds")
)

saveRDS(
  age_meta_final,
  file.path(output_dir, "data_processed/final_analysis/TCGA_BRCA_age_metadata_final.rds")
)

saveRDS(
  stage_meta_final,
  file.path(output_dir, "data_processed/final_analysis/TCGA_BRCA_stage_metadata_final.rds")
)

saveRDS(
  design_age_final,
  file.path(output_dir, "data_processed/final_analysis/TCGA_BRCA_age_design_final.rds")
)

saveRDS(
  design_stage_final,
  file.path(output_dir, "data_processed/final_analysis/TCGA_BRCA_stage_design_final.rds")
)
write.csv(
  data.frame(
    patient_id = final_qc_exclusions,
    exclusion_reason =
      "Concordant severe technical QC and multivariate expression-space abnormality"
  ),
  file.path(output_dir, "data_processed/final_analysis/TCGA_BRCA_final_QC_exclusions.csv"),
  row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "data_processed/final_analysis/sessionInfo_preprocessing_final.txt")
)
# ============================================================
# PHASE A — AGE ASSOCIATION
pipeline_event(ctx, "age_inference")
# A1. Final model validation
# ============================================================

summary(age_meta_final$age_years)

sd(age_meta_final$age_years)

anyNA(age_meta_final$age_years)

length(unique(age_meta_final$age_years))
hist(
  age_meta_final$age_years,
  breaks = 30,
  main = "Age at Diagnosis — Final Age Cohort",
  xlab = "Age at diagnosis (years)"
)
design_age_final
colnames(design_age_final)

head(design_age_final)

dim(design_age_final)

qr(design_age_final)$rank
identical(
  rownames(design_age_final),
  age_meta_final$patient_id
)

identical(
  colnames(voom_age_final$E),
  age_meta_final$patient_id
)

identical(
  colnames(voom_age_final$E),
  rownames(design_age_final)
)
rownames(design_age_final) <-
  age_meta_final$patient_id
identical(
  rownames(design_age_final),
  age_meta_final$patient_id
)

identical(
  colnames(voom_age_final$E),
  age_meta_final$patient_id
)

identical(
  colnames(voom_age_final$E),
  rownames(design_age_final)
)
head(
  data.frame(
    expression_sample =
      colnames(voom_age_final$E),
    
    metadata_sample =
      age_meta_final$patient_id,
    
    design_sample =
      rownames(design_age_final),
    
    age_metadata =
      age_meta_final$age_years,
    
    age_design =
      design_age_final[, "age_years"]
  ),
  10
)
identical(
  as.numeric(design_age_final[, "age_years"]),
  as.numeric(age_meta_final$age_years)
)
# ============================================================
# PHASE A — AGE ASSOCIATION
# A2. Gene-wise limma-voom regression
# ============================================================

fit_age <- lmFit(
  voom_age_final,
  design_age_final
)
pipeline_assert_fit(fit_age, voom_age_final, design_age_final)
class(fit_age)

dim(fit_age$coefficients)

colnames(fit_age$coefficients)

anyNA(fit_age$coefficients)
summary(
  fit_age$coefficients[, "age_years"]
)
range(
  fit_age$coefficients[, "age_years"]
)
# ============================================================
# A3. Empirical Bayes moderation
# ============================================================

fit_age_eb <- eBayes(fit_age)
class(fit_age_eb)

dim(fit_age_eb$coefficients)

anyNA(fit_age_eb$t)

anyNA(fit_age_eb$p.value)

summary(
  fit_age_eb$t[, "age_years"]
)

summary(
  fit_age_eb$p.value[, "age_years"]
)
range(
  fit_age_eb$p.value[, "age_years"]
)

sum(
  fit_age_eb$p.value[, "age_years"] < 0.05
)
# ============================================================
# A4. Age association results + BH-FDR
# ============================================================

age_results <- topTable(
  fit_age_eb,
  coef = "age_years",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)
pipeline_assert_results(age_results)
dim(age_results)

colnames(age_results)

head(age_results, 10)
sum(
  age_results$adj.P.Val < 0.05
)
sum(age_results$adj.P.Val < 0.01)

sum(age_results$adj.P.Val < 0.001)
summary(age_results$adj.P.Val)

range(age_results$adj.P.Val)
age_sig <- age_results[
  age_results$adj.P.Val < 0.05,
]

nrow(age_sig)

sum(age_sig$logFC > 0)

sum(age_sig$logFC < 0)
age_results$beta_age_per_year <-
  age_results$logFC

age_results$beta_age_per_10y <-
  age_results$logFC * 10
sum(age_sig$logFC > 0)

sum(age_sig$logFC < 0)
# ============================================================
# A4.1 — Construct interpretable Age results table
# ============================================================

age_results$ensembl_gene_id <- rownames(age_results)

age_results$beta_age_per_year <-
  age_results$logFC

age_results$beta_age_per_10y <-
  age_results$logFC * 10

age_results$direction <-
  ifelse(
    age_results$beta_age_per_year > 0,
    "Positive",
    "Negative"
  )

age_results$FDR_significant <-
  age_results$adj.P.Val < 0.05
age_results_final <- age_results[
  ,
  c(
    "ensembl_gene_id",
    "beta_age_per_year",
    "beta_age_per_10y",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B",
    "direction",
    "FDR_significant"
  )
]
dim(age_results_final)

head(age_results_final, 10)

table(age_results_final$FDR_significant)

table(
  age_results_final$direction[
    age_results_final$FDR_significant
  ]
)
age_results_final$ensembl_gene_id_base <-
  sub(
    "\\..*$",
    "",
    age_results_final$ensembl_gene_id
  )
head(
  age_results_final[
    ,
    c(
      "ensembl_gene_id",
      "ensembl_gene_id_base"
    )
  ],
  10
)
age_sig_final <-
  age_results_final[
    age_results_final$FDR_significant,
  ]

summary(
  abs(age_sig_final$beta_age_per_10y)
)
sum(
  abs(age_sig_final$beta_age_per_10y) >= 0.1
)

sum(
  abs(age_sig_final$beta_age_per_10y) >= 0.2
)

sum(
  abs(age_sig_final$beta_age_per_10y) >= 0.5
)
head(
  age_results_final[
    order(age_results_final$adj.P.Val),
  ],
  10
)
head(
  age_sig_final[
    order(
      abs(age_sig_final$beta_age_per_10y),
      decreasing = TRUE
    ),
  ],
  10
)
top_by_fdr <- head(
  age_results_final[
    order(age_results_final$adj.P.Val),
  ],
  10
)

top_by_effect <- head(
  age_sig_final[
    order(
      abs(age_sig_final$beta_age_per_10y),
      decreasing = TRUE
    ),
  ],
  10
)

age_candidates <- unique(
  c(
    top_by_fdr$ensembl_gene_id,
    top_by_effect$ensembl_gene_id
  )
)

length(age_candidates)
age_cor_results <- lapply(
  age_candidates,
  function(gene_id) {
    
    expr <- voom_age_final$E[
      gene_id,
    ]
    
    cor_test <- cor.test(
      age_meta_final$age_years,
      expr,
      method = "pearson"
    )
    
    data.frame(
      ensembl_gene_id = gene_id,
      correlation_r = unname(cor_test$estimate),
      correlation_p = cor_test$p.value
    )
  }
)

age_cor_results <- do.call(
  rbind,
  age_cor_results
)
age_cor_results
all(
  sign(age_cor_results$correlation_r) ==
    sign(
      age_results_final[
        match(
          age_cor_results$ensembl_gene_id,
          age_results_final$ensembl_gene_id
        ),
        "beta_age_per_year"
      ]
    )
)
ls()

grep(
  "gene|annot|feature",
  ls(),
  value = TRUE,
  ignore.case = TRUE
)
class(raw_count_matrix)
dim(raw_count_matrix)

head(
  rownames(raw_count_matrix)
)
age_cor_results <- merge(
  age_cor_results,
  age_results_final[
    ,
    c(
      "ensembl_gene_id",
      "beta_age_per_year",
      "beta_age_per_10y",
      "AveExpr",
      "adj.P.Val"
    )
  ],
  by = "ensembl_gene_id",
  all.x = TRUE
)

age_cor_results <- age_cor_results[
  order(abs(age_cor_results$correlation_r),
        decreasing = TRUE),
]

age_cor_results
range(abs(age_cor_results$correlation_r))

cor(
  age_cor_results$correlation_r,
  age_cor_results$beta_age_per_year,
  method = "pearson"
)
# ============================================================
# PHASE A — AGE ASSOCIATION
# A6. Inspect original GDC STAR gene annotation
# ============================================================

class(example_star)

dim(example_star)

colnames(example_star)

head(example_star, 10)
class(test_star)

dim(test_star)

colnames(test_star)

head(test_star, 10)
# ============================================================
# A6.1 — Build STAR-derived gene annotation table
# ============================================================

gene_annotation <- example_star[
  grepl("^ENSG", example_star$gene_id),
  c(
    "gene_id",
    "gene_name",
    "gene_type"
  )
]

dim(gene_annotation)

head(gene_annotation)

tail(gene_annotation)
length(unique(gene_annotation$gene_id))

anyDuplicated(gene_annotation$gene_id)

identical(
  gene_annotation$gene_id,
  rownames(raw_count_matrix)
)
# ============================================================
# A6.2 — Annotate Age association results
# ============================================================

age_results_annotated <- merge(
  age_results_final,
  gene_annotation,
  by.x = "ensembl_gene_id",
  by.y = "gene_id",
  all.x = TRUE,
  sort = FALSE
)
age_results_annotated <- age_results_annotated[
  match(
    age_results_final$ensembl_gene_id,
    age_results_annotated$ensembl_gene_id
  ),
]
dim(age_results_annotated)

identical(
  age_results_annotated$ensembl_gene_id,
  age_results_final$ensembl_gene_id
)

sum(is.na(age_results_annotated$gene_name))

sum(is.na(age_results_annotated$gene_type))

head(
  age_results_annotated[
    ,
    c(
      "ensembl_gene_id",
      "gene_name",
      "gene_type",
      "beta_age_per_10y",
      "adj.P.Val",
      "direction"
    )
  ],
  15
)
sort(
  table(age_results_annotated$gene_type),
  decreasing = TRUE
)
sort(
  table(
    age_results_annotated$gene_type[
      age_results_annotated$FDR_significant
    ]
  ),
  decreasing = TRUE
)
# ============================================================
# A6.3 — Data-driven candidate prioritization
# ============================================================

age_sig_annotated <- age_results_annotated[
  age_results_annotated$FDR_significant,
]

# Top statistical evidence
top_age_fdr <- head(
  age_sig_annotated[
    order(age_sig_annotated$adj.P.Val),
  ],
  20
)

# Top effect magnitude
top_age_effect <- head(
  age_sig_annotated[
    order(
      abs(age_sig_annotated$beta_age_per_10y),
      decreasing = TRUE
    ),
  ],
  20
)

# Protein-coding subset
age_sig_protein <- age_sig_annotated[
  age_sig_annotated$gene_type == "protein_coding",
]

top_age_protein_fdr <- head(
  age_sig_protein[
    order(age_sig_protein$adj.P.Val),
  ],
  20
)

top_age_protein_effect <- head(
  age_sig_protein[
    order(
      abs(age_sig_protein$beta_age_per_10y),
      decreasing = TRUE
    ),
  ],
  20
)
top_age_fdr[
  ,
  c(
    "gene_name",
    "gene_type",
    "beta_age_per_10y",
    "AveExpr",
    "adj.P.Val",
    "direction"
  )
]
top_age_effect[
  ,
  c(
    "gene_name",
    "gene_type",
    "beta_age_per_10y",
    "AveExpr",
    "adj.P.Val",
    "direction"
  )
]
top_age_protein_fdr[
  ,
  c(
    "gene_name",
    "beta_age_per_10y",
    "AveExpr",
    "adj.P.Val",
    "direction"
  )
]

top_age_protein_effect[
  ,
  c(
    "gene_name",
    "beta_age_per_10y",
    "AveExpr",
    "adj.P.Val",
    "direction"
  )
]
gene_type_summary <- data.frame(
  tested = table(age_results_annotated$gene_type),
  significant = table(
    factor(
      age_results_annotated$gene_type[
        age_results_annotated$FDR_significant
      ],
      levels = names(table(age_results_annotated$gene_type))
    )
  )
)

gene_type_summary$significant_fraction <-
  gene_type_summary$significant.Freq /
  gene_type_summary$tested.Freq

gene_type_summary <- gene_type_summary[
  order(
    gene_type_summary$significant_fraction,
    decreasing = TRUE
  ),
]

head(gene_type_summary, 20)
# ============================================================
# A7 — Publication-quality visualization of Age associations
# ============================================================

library(ggplot2)

age_figures_dir <- file.path(
  output_dir,
  "figures",
  "age_analysis"
)

dir.create(
  age_figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
theme_age_project <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    plot.caption = element_text(
      size = 9,
      hjust = 0
    )
  )
#Prepare plotting data
age_volcano <- age_results_annotated

age_volcano$neg_log10_FDR <- -log10(
  age_volcano$adj.P.Val
)

age_volcano$association_class <- "Not significant"

age_volcano$association_class[
  age_volcano$adj.P.Val < 0.05 &
    age_volcano$beta_age_per_10y > 0
] <- "Positive"

age_volcano$association_class[
  age_volcano$adj.P.Val < 0.05 &
    age_volcano$beta_age_per_10y < 0
] <- "Negative"

table(age_volcano$association_class)
stopifnot(
  sum(age_volcano$association_class == "Positive") == 4367,
  sum(age_volcano$association_class == "Negative") == 4696
)
volcano_label_genes <- c(
  "ESR1",
  "SOX10",
  "GABRP",
  "KRT5",
  "KRT17",
  "DSC2",
  "TNK1",
  "CLEC3A"
)

age_volcano_labels <- age_volcano[
  age_volcano$gene_name %in% volcano_label_genes,
]

age_volcano_labels[
  ,
  c(
    "gene_name",
    "beta_age_per_10y",
    "AveExpr",
    "adj.P.Val"
  )
]
library(ggrepel)
# ============================================================
# A7.1 — Final transcriptome-wide Age association plot
# ============================================================

fig_age_volcano <- ggplot(
  age_volcano,
  aes(
    x = beta_age_per_10y,
    y = neg_log10_FDR
  )
) +
  
  # All genes
  geom_point(
    aes(color = association_class),
    alpha = 0.42,
    size = 1.15
  ) +
  
  # BH-FDR threshold
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  # Zero effect
  geom_vline(
    xintercept = 0,
    linewidth = 0.45
  ) +
  
  # Selected biologically interpretable candidates
  geom_text_repel(
    data = age_volcano_labels,
    aes(label = gene_name),
    size = 3.5,
    fontface = "plain",
    box.padding = 0.55,
    point.padding = 0.35,
    min.segment.length = 0,
    segment.size = 0.3,
    max.overlaps = Inf,
    seed = 42,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      "Negative" = "#3B6FB6",
      "Not significant" = "grey78",
      "Positive" = "#D6534D"
    ),
    breaks = c(
      "Negative",
      "Not significant",
      "Positive"
    )
  ) +
  
  scale_x_continuous(
    breaks = seq(-0.8, 0.8, by = 0.2),
    expand = expansion(mult = c(0.08, 0.10))
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  labs(
    title = "Transcriptome-wide association with age at diagnosis",
    subtitle = "TCGA-BRCA primary tumors (n = 1,072; 25,551 genes tested)",
    x = expression(
      beta[Age]~"per 10 years (voom log"[2]*"-CPM)"
    ),
    y = expression(
      -log[10]*"(BH-FDR)"
    ),
    color = NULL,
    caption = "Dashed line: BH-FDR = 0.05."
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      margin = margin(b = 5)
    ),
    
    plot.subtitle = element_text(
      size = 11,
      margin = margin(b = 10)
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    legend.position = "top",
    
    legend.justification = "left",
    
    legend.margin = margin(
      t = 0,
      r = 0,
      b = 5,
      l = 0
    ),
    
    legend.box.spacing = unit(
      2,
      "pt"
    ),
    
    plot.caption = element_text(
      size = 9,
      hjust = 0,
      margin = margin(t = 8)
    ),
    
    plot.margin = margin(
      t = 12,
      r = 25,
      b = 12,
      l = 12
    )
  )

fig_age_volcano
ggsave(
  filename = file.path(
    age_figures_dir,
    "Figure_Age_01_transcriptome_wide_association.png"
  ),
  plot = fig_age_volcano,
  width = 9,
  height = 6.8,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    age_figures_dir,
    "Figure_Age_01_transcriptome_wide_association.pdf"
  ),
  plot = fig_age_volcano,
  width = 9,
  height = 6.8,
  units = "in",
  bg = "white"
)

file.exists(
  file.path(
    age_figures_dir,
    "Figure_Age_01_transcriptome_wide_association.png"
  )
)

file.info(
  file.path(
    age_figures_dir,
    "Figure_Age_01_transcriptome_wide_association.png"
  )
)$size
# ============================================================
# A7.2 — Selected gene Age-expression scatterplots
# ============================================================

scatter_genes <- c(
  "ESR1",
  "SOX10",
  "GABRP",
  "KRT5"
)

scatter_gene_info <- age_results_annotated[
  age_results_annotated$gene_name %in% scatter_genes,
  c(
    "ensembl_gene_id",
    "gene_name",
    "beta_age_per_10y",
    "AveExpr",
    "adj.P.Val"
  )
]

# Put genes in our intended display order
scatter_gene_info <- scatter_gene_info[
  match(scatter_genes, scatter_gene_info$gene_name),
]

scatter_gene_info

stopifnot(
  nrow(scatter_gene_info) == 4,
  !any(is.na(scatter_gene_info$ensembl_gene_id))
)
scatter_expression <- voom_age_final$E[
  scatter_gene_info$ensembl_gene_id,
  ,
  drop = FALSE
]

dim(scatter_expression)

identical(
  colnames(scatter_expression),
  age_meta_final$patient_id
)
scatter_data <- do.call(
  rbind,
  lapply(
    seq_len(nrow(scatter_gene_info)),
    function(i) {
      
      data.frame(
        patient_id = age_meta_final$patient_id,
        age_years = age_meta_final$age_years,
        gene_name = scatter_gene_info$gene_name[i],
        expression = as.numeric(
          scatter_expression[i, ]
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(scatter_data) <- NULL

dim(scatter_data)

head(scatter_data)
table(scatter_data$gene_name)

sum(is.na(scatter_data$age_years))
sum(is.na(scatter_data$expression))
pearson_summary <- do.call(
  rbind,
  lapply(
    scatter_genes,
    function(g) {
      
      tmp <- scatter_data[
        scatter_data$gene_name == g,
        ,
        drop = FALSE
      ]
      
      data.frame(
        gene_name = g,
        pearson_r = cor(
          tmp$age_years,
          tmp$expression,
          method = "pearson"
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

pearson_summary
scatter_stats <- merge(
  scatter_gene_info,
  pearson_summary,
  by = "gene_name",
  sort = FALSE
)

scatter_stats <- scatter_stats[
  match(scatter_genes, scatter_stats$gene_name),
]

scatter_stats
all(
  sign(scatter_stats$beta_age_per_10y) ==
    sign(scatter_stats$pearson_r)
)
scatter_stats$annotation <- sprintf(
  "\u03B210y = %+.2f\nFDR = %.2e\nPearson r = %+.2f",
  scatter_stats$beta_age_per_10y,
  scatter_stats$adj.P.Val,
  scatter_stats$pearson_r
)

scatter_stats[
  ,
  c(
    "gene_name",
    "beta_age_per_10y",
    "adj.P.Val",
    "pearson_r",
    "annotation"
  )
]
scatter_data$gene_name <- factor(
  scatter_data$gene_name,
  levels = scatter_genes
)

scatter_stats$gene_name <- factor(
  scatter_stats$gene_name,
  levels = scatter_genes
)
# A7.2 — CANONICAL refined scatterplot; final labels and one PNG/PDF export.
scatter_stats$annotation_short <- sprintf(
  "\u03B210y = %+.2f\nFDR = %.2e\nr = %+.2f",
  scatter_stats$beta_age_per_10y,
  scatter_stats$adj.P.Val,
  scatter_stats$pearson_r
)

fig_age_scatter_final <- ggplot(
  scatter_data,
  aes(x = age_years, y = expression)
) +
  
  geom_point(
    alpha = 0.18,
    size = 0.8
  ) +
  
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.8
  ) +
  
  facet_wrap(
    ~ gene_name,
    scales = "free_y",
    ncol = 2
  ) +
  
  geom_label(
    data = scatter_stats,
    aes(
      x = -Inf,
      y = Inf,
      label = annotation_short
    ),
    inherit.aes = FALSE,
    hjust = -0.08,
    vjust = 1.10,
    size = 2.7,
    lineheight = 0.95,
    linewidth = 0.15,
    label.padding = unit(0.10, "lines"),
    fill = "white",
    alpha = 0.78
  ) +
  
  labs(
    title = "Selected gene-expression associations with age at diagnosis",
    subtitle = "TCGA-BRCA primary tumors (n = 1,072)",
    x = "Age at diagnosis (years)",
    y = expression("Voom log"[2]*"-CPM expression"),
    caption = paste0(
      "\u03B210y and FDR: primary limma-voom model; ",
      "r: descriptive Pearson correlation; ",
      "lines and 95% CIs: descriptive ordinary linear trends."
    )
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      size = 11,
      margin = margin(b = 10)
    ),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    panel.spacing = unit(1.3, "lines"),
    plot.caption = element_text(
      size = 8.3,
      hjust = 0,
      margin = margin(t = 9)
    ),
    plot.margin = margin(
      12, 22, 12, 12
    )
  )

fig_age_scatter_final
ggsave(
  file.path(
    age_figures_dir,
    "Figure_Age_02_selected_gene_scatterplots.png"
  ),
  fig_age_scatter_final,
  width = 9,
  height = 7.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(
    age_figures_dir,
    "Figure_Age_02_selected_gene_scatterplots.pdf"
  ),
  fig_age_scatter_final,
  width = 9,
  height = 7.5,
  units = "in",
  bg = "white"
)
# ============================================================
# A7.3 — Ranked Age effect sizes
# Protein-coding, FDR-significant genes
# ============================================================

age_effect_candidates <- age_results_annotated[
  age_results_annotated$gene_type == "protein_coding" &
    age_results_annotated$FDR_significant,
  ,
  drop = FALSE
]

dim(age_effect_candidates)

table(age_effect_candidates$direction)
age_effect_candidates$abs_beta_age_per_10y <-
  abs(age_effect_candidates$beta_age_per_10y)

age_effect_candidates <- age_effect_candidates[
  order(
    age_effect_candidates$abs_beta_age_per_10y,
    decreasing = TRUE
  ),
]

head(
  age_effect_candidates[
    ,
    c(
      "gene_name",
      "beta_age_per_10y",
      "AveExpr",
      "adj.P.Val",
      "direction"
    )
  ],
  20
)
top_n_effect <- 30

age_effect_top30 <- head(
  age_effect_candidates,
  top_n_effect
)

dim(age_effect_top30)

age_effect_top30[
  ,
  c(
    "gene_name",
    "beta_age_per_10y",
    "AveExpr",
    "adj.P.Val",
    "direction"
  )
]
stopifnot(
  nrow(age_effect_top30) == 30,
  all(age_effect_top30$gene_type == "protein_coding"),
  all(age_effect_top30$adj.P.Val < 0.05)
)

range(age_effect_top30$beta_age_per_10y)

table(age_effect_top30$direction)
age_effect_top30$gene_name_plot <- factor(
  age_effect_top30$gene_name,
  levels = age_effect_top30$gene_name[
    order(age_effect_top30$beta_age_per_10y)
  ]
)
age_effect_top30[
  order(age_effect_top30$beta_age_per_10y),
  c(
    "gene_name",
    "beta_age_per_10y",
    "direction"
  )
]
# ============================================================
# Figure 3 — Ranked Age effect sizes
# ============================================================

fig_age_effect_ranking <- ggplot(
  age_effect_top30,
  aes(
    x = beta_age_per_10y,
    y = gene_name_plot
  )
) +
  
  geom_segment(
    aes(
      x = 0,
      xend = beta_age_per_10y,
      y = gene_name_plot,
      yend = gene_name_plot,
      color = direction
    ),
    linewidth = 0.7
  ) +
  
  geom_point(
    aes(color = direction),
    size = 2.6
  ) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.45
  ) +
  
  scale_color_manual(
    values = c(
      "Negative" = "#3B6FB6",
      "Positive" = "#D6534D"
    )
  ) +
  
  labs(
    title = "Largest estimated Age-associated expression effects",
    subtitle = paste0(
      "Top 30 FDR-significant protein-coding genes ",
      "ranked by absolute effect size"
    ),
    x = expression(
      beta[Age]~"per 10 years (voom log"[2]*"-CPM)"
    ),
    y = NULL,
    color = NULL,
    caption = paste0(
      "Genes: protein-coding and BH-FDR < 0.05; ",
      "ranking based on absolute \u03B2Age per 10 years."
    )
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      margin = margin(b = 5)
    ),
    
    plot.subtitle = element_text(
      size = 10.5,
      margin = margin(b = 10)
    ),
    
    axis.title.x = element_text(
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 9
    ),
    
    legend.position = "top",
    
    plot.caption = element_text(
      size = 8.5,
      hjust = 0,
      margin = margin(t = 8)
    ),
    
    plot.margin = margin(
      t = 12,
      r = 20,
      b = 12,
      l = 12
    )
  )

fig_age_effect_ranking
# ============================================================
# A7.3 — Final validation of ranked candidates
# ============================================================

age_effect_top30[
  ,
  c(
    "gene_name",
    "ensembl_gene_id",
    "beta_age_per_10y",
    "AveExpr",
    "adj.P.Val",
    "direction"
  )
]

# Number in each direction
table(age_effect_top30$direction)

# Check gene symbols
anyDuplicated(age_effect_top30$gene_name)
sum(is.na(age_effect_top30$gene_name))
sum(age_effect_top30$gene_name == "")

# Expression-abundance range
summary(age_effect_top30$AveExpr)

# Which top-effect genes have relatively low average expression?
age_effect_top30[
  order(age_effect_top30$AveExpr),
  c(
    "gene_name",
    "beta_age_per_10y",
    "AveExpr",
    "adj.P.Val"
  )
]
fig_age_effect_ranking_final <- fig_age_effect_ranking +
  labs(
    title = "Largest estimated age-associated expression effects",
    caption = paste0(
      "Eligible genes: protein-coding and BH-FDR < 0.05; ",
      "top 30 ranked by absolute \u03B2Age per 10 years."
    )
  )

ggsave(
  file.path(
    age_figures_dir,
    "Figure_Age_03_ranked_effect_sizes.png"
  ),
  fig_age_effect_ranking_final,
  width = 8.5,
  height = 7.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(
    age_figures_dir,
    "Figure_Age_03_ranked_effect_sizes.pdf"
  ),
  fig_age_effect_ranking_final,
  width = 8.5,
  height = 7.5,
  units = "in",
  bg = "white"
)
#A8 — Preliminary Biological Interpretation of Age Associations
# ============================================================
# A8.1 — Freeze primary Age findings for interpretation
# ============================================================

age_interpretation_summary <- data.frame(
  metric = c(
    "Patients",
    "Genes tested",
    "FDR-significant genes",
    "FDR-significant positive",
    "FDR-significant negative",
    "FDR-significant protein-coding",
    "Top-30 effect genes: positive",
    "Top-30 effect genes: negative"
  ),
  value = c(
    nrow(age_meta_final),
    nrow(age_results_annotated),
    sum(age_results_annotated$adj.P.Val < 0.05),
    sum(
      age_results_annotated$adj.P.Val < 0.05 &
        age_results_annotated$direction == "Positive"
    ),
    sum(
      age_results_annotated$adj.P.Val < 0.05 &
        age_results_annotated$direction == "Negative"
    ),
    sum(
      age_results_annotated$adj.P.Val < 0.05 &
        age_results_annotated$gene_type == "protein_coding"
    ),
    sum(age_effect_top30$direction == "Positive"),
    sum(age_effect_top30$direction == "Negative")
  )
)

age_interpretation_summary
# ============================================================
# A8.2 — Candidate interpretation set
# ============================================================

age_interpretation_genes <- c(
  # Positive, high-interest
  "ESR1",
  "PIP",
  "CYP4Z1",
  "TNK1",
  
  # Negative, basal/epithelial-associated candidates
  "SOX10",
  "GABRP",
  "KRT5",
  "KRT14",
  "KRT17",
  "DSC3",
  "DSG3",
  
  # Other strong negative associations
  "SFRP1",
  "SLC34A2",
  "SERPINB5",
  
  # Extreme effect but abundance caution
  "CLEC3A",
  "DCX",
  "DLK1"
)

age_interpretation_table <- age_results_annotated[
  age_results_annotated$gene_name %in% age_interpretation_genes,
  c(
    "gene_name",
    "gene_type",
    "beta_age_per_10y",
    "AveExpr",
    "adj.P.Val",
    "direction"
  )
]

age_interpretation_table <- age_interpretation_table[
  match(
    age_interpretation_genes,
    age_interpretation_table$gene_name
  ),
]

age_interpretation_table
# ============================================================
# A8.3.1 — Contextualize the ESR1-associated signal
# ============================================================

luminal_context_genes <- c(
  "ESR1",
  "PGR",
  "GATA3",
  "FOXA1",
  "XBP1"
)

luminal_context <- age_results_annotated[
  age_results_annotated$gene_name %in% luminal_context_genes,
  c(
    "gene_name",
    "beta_age_per_10y",
    "AveExpr",
    "adj.P.Val",
    "direction"
  )
]

luminal_context <- luminal_context[
  match(
    luminal_context_genes,
    luminal_context$gene_name
  ),
]

luminal_context
# Stage inference reuses canonical STEP 8 objects; no second cohort/filter/TMM/voom.
stopifnot(
  identical(stage_meta_final$patient_id, colnames(raw_count_stage_final)),
  identical(stage_meta_final$patient_id, rownames(design_stage_final)),
  identical(stage_meta_final$patient_id, colnames(voom_stage_final$E)),
  !any(final_qc_exclusions %in% stage_meta_final$patient_id),
  !anyNA(stage_meta_final$stage_group),
  qr(design_stage_final)$rank == ncol(design_stage_final)
)
effective_lib_sizes_stage <-
  dge_stage_final_norm$samples$lib.size * dge_stage_final_norm$samples$norm.factors
# ============================================================
# S1 — FINAL STAGE DATASET & DESIGN FINGERPRINT
# ============================================================

stage_fingerprint <- c(
  patients = nrow(stage_meta_final),
  genes_tested = nrow(voom_stage_final$E),
  stage_I = sum(stage_meta_final$stage_group == "Stage I"),
  stage_II = sum(stage_meta_final$stage_group == "Stage II"),
  stage_III = sum(stage_meta_final$stage_group == "Stage III"),
  stage_IV = sum(stage_meta_final$stage_group == "Stage IV"),
  missing_stage = sum(is.na(stage_meta_final$stage_group)),
  design_columns = ncol(design_stage_final),
  design_rank = qr(design_stage_final)$rank
)

print(stage_fingerprint)

stopifnot(
  stage_fingerprint["patients"] == 1063,
  stage_fingerprint["genes_tested"] == 31157,
  stage_fingerprint["stage_I"] == 179,
  stage_fingerprint["stage_II"] == 615,
  stage_fingerprint["stage_III"] == 249,
  stage_fingerprint["stage_IV"] == 20,
  stage_fingerprint["missing_stage"] == 0,
  stage_fingerprint["design_columns"] == 4,
  stage_fingerprint["design_rank"] == 4,
  
  identical(
    colnames(voom_stage_final$E),
    stage_meta_final$patient_id
  ),
  
  identical(
    colnames(voom_stage_final$E),
    rownames(design_stage_final)
  )
)

cat("\nS1 VALIDATION PASSED — STAGE DATASET FROZEN\n")
# ============================================================
pipeline_event(ctx, "stage_inference")
# S2 — Fit primary Stage gene-wise linear models
# ============================================================

fit_stage <- lmFit(
  voom_stage_final,
  design_stage_final
)
pipeline_assert_fit(fit_stage, voom_stage_final, design_stage_final)

# Basic structure
class(fit_stage)

dim(fit_stage$coefficients)

colnames(fit_stage$coefficients)

head(fit_stage$coefficients)
# ============================================================
# S2.2 — Validate fitted model
# ============================================================

# Dimensions
dim(fit_stage$coefficients)

# Missing coefficients
sum(is.na(fit_stage$coefficients))

# Non-finite coefficients
sum(!is.finite(fit_stage$coefficients))

# Residual degrees of freedom
summary(fit_stage$df.residual)

# Residual standard deviations
summary(fit_stage$sigma)

# Check gene identity/order
identical(
  rownames(fit_stage$coefficients),
  rownames(voom_stage_final$E)
)
# ============================================================
# S2.3 — Inspect Stage coefficient distributions
# ============================================================

stage_coef_summary <- data.frame(
  contrast = c(
    "Stage II - Stage I",
    "Stage III - Stage I",
    "Stage IV - Stage I"
  ),
  
  min = c(
    min(fit_stage$coefficients[, "stage_groupStage II"]),
    min(fit_stage$coefficients[, "stage_groupStage III"]),
    min(fit_stage$coefficients[, "stage_groupStage IV"])
  ),
  
  median = c(
    median(fit_stage$coefficients[, "stage_groupStage II"]),
    median(fit_stage$coefficients[, "stage_groupStage III"]),
    median(fit_stage$coefficients[, "stage_groupStage IV"])
  ),
  
  max = c(
    max(fit_stage$coefficients[, "stage_groupStage II"]),
    max(fit_stage$coefficients[, "stage_groupStage III"]),
    max(fit_stage$coefficients[, "stage_groupStage IV"])
  )
)

stage_coef_summary
# ============================================================
# S3.1 — Empirical Bayes moderation
# ============================================================

fit_stage_eb <- eBayes(fit_stage)

class(fit_stage_eb)

dim(fit_stage_eb$t)

sum(is.na(fit_stage_eb$t))
sum(!is.finite(fit_stage_eb$t))

summary(fit_stage_eb$df.total)

summary(fit_stage_eb$s2.post)
stage_coef_names <- c(
  "stage_groupStage II",
  "stage_groupStage III",
  "stage_groupStage IV"
)

stage_coef_names
# ============================================================
# S3.2 — Omnibus moderated F-test for Stage
# ============================================================

stage_omnibus <- topTable(
  fit_stage_eb,
  coef = stage_coef_names,
  number = Inf,
  adjust.method = "BH",
  sort.by = "F"
)
pipeline_assert_results(stage_omnibus)

dim(stage_omnibus)

names(stage_omnibus)

head(stage_omnibus, 10)
# ============================================================
# S3.3 — Validate omnibus statistics
# ============================================================

sum(is.na(stage_omnibus$F))
sum(is.na(stage_omnibus$P.Value))
sum(is.na(stage_omnibus$adj.P.Val))

range(stage_omnibus$F)
range(stage_omnibus$P.Value)
range(stage_omnibus$adj.P.Val)

summary(stage_omnibus$F)

# topTable returns syntactically valid names; fit/design names retain spaces.
stopifnot(all(make.names(stage_coef_names) %in% names(stage_omnibus)))
head(
  stage_omnibus[
    ,
    c(
      make.names(stage_coef_names),
      "AveExpr",
      "F",
      "P.Value",
      "adj.P.Val"
    )
  ],
  10
)
identical(
  sort(rownames(stage_omnibus)),
  sort(rownames(voom_stage_final$E))
)
head(
  stage_omnibus[
    ,
    c(
      "stage_groupStage.II",
      "stage_groupStage.III",
      "stage_groupStage.IV",
      "AveExpr",
      "F",
      "P.Value",
      "adj.P.Val"
    )
  ],
  10
)
# ============================================================
# S4.1 — Define Stage-associated genes by BH-FDR
# ============================================================

stage_omnibus$ensembl_gene_id <- rownames(stage_omnibus)

stage_omnibus$FDR_significant <-
  stage_omnibus$adj.P.Val < 0.05

table(stage_omnibus$FDR_significant)

sum(stage_omnibus$FDR_significant)

mean(stage_omnibus$FDR_significant) * 100
# ============================================================
# S4.2 — Significance summary
# ============================================================

stage_fdr_summary <- c(
  FDR_0.05 = sum(stage_omnibus$adj.P.Val < 0.05),
  FDR_0.01 = sum(stage_omnibus$adj.P.Val < 0.01),
  FDR_0.001 = sum(stage_omnibus$adj.P.Val < 0.001)
)

stage_fdr_summary
# ============================================================
# S4.3 — Inspect strongest omnibus Stage associations
# ============================================================

head(
  stage_omnibus[
    order(stage_omnibus$adj.P.Val),
    c(
      "ensembl_gene_id",
      "stage_groupStage.II",
      "stage_groupStage.III",
      "stage_groupStage.IV",
      "AveExpr",
      "F",
      "P.Value",
      "adj.P.Val",
      "FDR_significant"
    )
  ],
  20
)
# ============================================================
# S5.1 — CANONICAL six pairwise Stage contrasts (explicit matrix)
# ============================================================
stopifnot(identical(
  colnames(design_stage_final),
  c("(Intercept)", "stage_groupStage II", "stage_groupStage III", "stage_groupStage IV")
))

stage_contrasts <- matrix(
  c(
    # II-I  III-I  IV-I  III-II IV-II IV-III
    0,     0,     0,      0,     0,     0,   # Intercept
    1,     0,     0,     -1,    -1,     0,   # Stage II
    0,     1,     0,      1,     0,    -1,   # Stage III
    0,     0,     1,      0,     1,     1    # Stage IV
  ),
  nrow = 4,
  byrow = TRUE
)

rownames(stage_contrasts) <- colnames(design_stage_final)

colnames(stage_contrasts) <- c(
  "II_vs_I",
  "III_vs_I",
  "IV_vs_I",
  "III_vs_II",
  "IV_vs_II",
  "IV_vs_III"
)

stage_contrasts
# ============================================================
# S5.2 — Fit all pairwise Stage contrasts
# ============================================================

fit_stage_contrasts <- contrasts.fit(
  fit_stage,
  stage_contrasts
)

fit_stage_contrasts_eb <- eBayes(
  fit_stage_contrasts
)

# Basic validation
class(fit_stage_contrasts_eb)

dim(fit_stage_contrasts_eb$coefficients)

colnames(fit_stage_contrasts_eb$coefficients)

sum(is.na(fit_stage_contrasts_eb$coefficients))
sum(!is.finite(fit_stage_contrasts_eb$coefficients))

sum(is.na(fit_stage_contrasts_eb$t))
sum(!is.finite(fit_stage_contrasts_eb$t))

summary(fit_stage_contrasts_eb$df.total)
# ============================================================
# S5.3 — Verify reference contrasts against original model
# ============================================================

all.equal(
  fit_stage_contrasts_eb$coefficients[, "II_vs_I"],
  fit_stage$coefficients[, "stage_groupStage II"]
)

all.equal(
  fit_stage_contrasts_eb$coefficients[, "III_vs_I"],
  fit_stage$coefficients[, "stage_groupStage III"]
)

all.equal(
  fit_stage_contrasts_eb$coefficients[, "IV_vs_I"],
  fit_stage$coefficients[, "stage_groupStage IV"]
)
all.equal(
  fit_stage_contrasts_eb$coefficients[, "III_vs_II"],
  fit_stage$coefficients[, "stage_groupStage III"] -
    fit_stage$coefficients[, "stage_groupStage II"]
)

all.equal(
  fit_stage_contrasts_eb$coefficients[, "IV_vs_III"],
  fit_stage$coefficients[, "stage_groupStage IV"] -
    fit_stage$coefficients[, "stage_groupStage III"]
)
# ============================================================
# S5.4 — Extract pairwise contrast results
# ============================================================

stage_contrast_results <- lapply(
  colnames(stage_contrasts),
  function(x) {
    
    res <- topTable(
      fit_stage_contrasts_eb,
      coef = x,
      number = Inf,
      adjust.method = "BH",
      sort.by = "P"
    )
    
    res$ensembl_gene_id <- rownames(res)
    res$contrast <- x
    
    res
  }
)

names(stage_contrast_results) <- colnames(stage_contrasts)

sapply(stage_contrast_results, dim)
# ============================================================
invisible(lapply(stage_contrast_results, pipeline_assert_results))
# S5.5 — Pairwise contrast summary
# ============================================================

stage_contrast_summary <- do.call(
  rbind,
  lapply(
    names(stage_contrast_results),
    function(x) {
      
      res <- stage_contrast_results[[x]]
      
      data.frame(
        contrast = x,
        FDR_0.05 = sum(res$adj.P.Val < 0.05),
        FDR_0.01 = sum(res$adj.P.Val < 0.01),
        positive_FDR_0.05 = sum(
          res$adj.P.Val < 0.05 & res$logFC > 0
        ),
        negative_FDR_0.05 = sum(
          res$adj.P.Val < 0.05 & res$logFC < 0
        ),
        min_logFC = min(res$logFC),
        max_logFC = max(res$logFC)
      )
    }
  )
)

rownames(stage_contrast_summary) <- NULL

stage_contrast_summary
# ============================================================
# S5.6 — Connect omnibus and pairwise results
# ============================================================

omnibus_sig_genes <- stage_omnibus$ensembl_gene_id[
  stage_omnibus$adj.P.Val < 0.05
]

pairwise_sig_sets <- lapply(
  stage_contrast_results,
  function(res) {
    res$ensembl_gene_id[res$adj.P.Val < 0.05]
  }
)

pairwise_union <- unique(
  unlist(pairwise_sig_sets)
)

# Omnibus genes also significant in >=1 pairwise contrast
omnibus_with_pairwise <- intersect(
  omnibus_sig_genes,
  pairwise_union
)

# Omnibus genes without any individually significant pairwise contrast
omnibus_only <- setdiff(
  omnibus_sig_genes,
  pairwise_union
)

c(
  omnibus_significant = length(omnibus_sig_genes),
  pairwise_union = length(pairwise_union),
  omnibus_with_pairwise = length(omnibus_with_pairwise),
  omnibus_only = length(omnibus_only)
)

length(omnibus_with_pairwise) /
  length(omnibus_sig_genes) * 100
# Are any pairwise-significant genes NOT significant
# in the primary omnibus test?

pairwise_not_omnibus <- setdiff(
  pairwise_union,
  omnibus_sig_genes
)

length(pairwise_not_omnibus)
# ============================================================
# S6.1 — Validate canonical GENCODE v36 annotation
# ============================================================

dim(gene_annotation)
names(gene_annotation)

head(gene_annotation)

sum(duplicated(gene_annotation$gene_id))
sum(is.na(gene_annotation$gene_name))
sum(is.na(gene_annotation$gene_type))

identical(
  gene_annotation$gene_id,
  rownames(raw_count_matrix)
)
# ============================================================
# S6.2 — Annotate omnibus Stage results
# ============================================================

annotation_idx <- match(
  stage_omnibus$ensembl_gene_id,
  gene_annotation$gene_id
)

sum(is.na(annotation_idx))

stage_results_annotated <- cbind(
  stage_omnibus,
  gene_name = gene_annotation$gene_name[annotation_idx],
  gene_type = gene_annotation$gene_type[annotation_idx]
)

dim(stage_results_annotated)

sum(is.na(stage_results_annotated$gene_name))
sum(is.na(stage_results_annotated$gene_type))

head(
  stage_results_annotated[
    ,
    c(
      "ensembl_gene_id",
      "gene_name",
      "gene_type",
      "stage_groupStage.II",
      "stage_groupStage.III",
      "stage_groupStage.IV",
      "AveExpr",
      "F",
      "adj.P.Val"
    )
  ],
  20
)
# ============================================================
# S6.3 — Biotype composition
# ============================================================

stage_sig <- stage_results_annotated[
  stage_results_annotated$FDR_significant,
]

dim(stage_sig)

sort(
  table(stage_sig$gene_type),
  decreasing = TRUE
)

stage_biotype_summary <- data.frame(
  gene_type = names(table(stage_results_annotated$gene_type)),
  tested = as.integer(table(stage_results_annotated$gene_type))
)

sig_counts <- table(stage_sig$gene_type)

stage_biotype_summary$significant <- as.integer(
  sig_counts[
    match(
      stage_biotype_summary$gene_type,
      names(sig_counts)
    )
  ]
)

stage_biotype_summary$significant[
  is.na(stage_biotype_summary$significant)
] <- 0

stage_biotype_summary$percent_significant <-
  100 *
  stage_biotype_summary$significant /
  stage_biotype_summary$tested

stage_biotype_summary <- stage_biotype_summary[
  order(
    stage_biotype_summary$significant,
    decreasing = TRUE
  ),
]

head(stage_biotype_summary, 20)
# ============================================================
# S6.4 — Top primary Stage-associated genes
# ============================================================

top_stage_omnibus <- stage_results_annotated[
  order(stage_results_annotated$adj.P.Val),
  c(
    "ensembl_gene_id",
    "gene_name",
    "gene_type",
    "stage_groupStage.II",
    "stage_groupStage.III",
    "stage_groupStage.IV",
    "AveExpr",
    "F",
    "P.Value",
    "adj.P.Val"
  )
]

head(top_stage_omnibus, 30)
# ============================================================
# S6.5 — Effect-size characterization of omnibus-significant genes
# ============================================================

stage_sig$max_abs_stage_effect <- apply(
  abs(
    stage_sig[
      ,
      c(
        "stage_groupStage.II",
        "stage_groupStage.III",
        "stage_groupStage.IV"
      )
    ]
  ),
  1,
  max
)

summary(stage_sig$max_abs_stage_effect)

# Protein-coding Stage-associated genes
stage_sig_pc <- stage_sig[
  stage_sig$gene_type == "protein_coding",
]

# Rank by maximum absolute Stage-vs-I coefficient
stage_sig_pc_effect_ranked <- stage_sig_pc[
  order(
    stage_sig_pc$max_abs_stage_effect,
    decreasing = TRUE
  ),
]

head(
  stage_sig_pc_effect_ranked[
    ,
    c(
      "ensembl_gene_id",
      "gene_name",
      "stage_groupStage.II",
      "stage_groupStage.III",
      "stage_groupStage.IV",
      "max_abs_stage_effect",
      "AveExpr",
      "F",
      "adj.P.Val"
    )
  ],
  30
)
stage_sig_pc_abundant <- stage_sig_pc[
  stage_sig_pc$AveExpr >= 1,
]

stage_sig_pc_abundant <- stage_sig_pc_abundant[
  order(
    stage_sig_pc_abundant$max_abs_stage_effect,
    decreasing = TRUE
  ),
]

head(
  stage_sig_pc_abundant[
    ,
    c(
      "gene_name",
      "stage_groupStage.II",
      "stage_groupStage.III",
      "stage_groupStage.IV",
      "max_abs_stage_effect",
      "AveExpr",
      "adj.P.Val"
    )
  ],
  30
)
# ============================================================
# S7.1 — Prepare omnibus Stage association overview
# ============================================================

stage_plot_data <- stage_results_annotated

stage_plot_data$max_abs_stage_effect <- apply(
  abs(
    stage_plot_data[
      ,
      c(
        "stage_groupStage.II",
        "stage_groupStage.III",
        "stage_groupStage.IV"
      )
    ]
  ),
  1,
  max
)

stage_plot_data$neg_log10_FDR <-
  -log10(stage_plot_data$adj.P.Val)

stage_plot_data$association_class <-
  ifelse(
    stage_plot_data$adj.P.Val < 0.05,
    "Stage-associated",
    "Not significant"
  )


table(stage_plot_data$association_class)

summary(stage_plot_data$max_abs_stage_effect)
summary(stage_plot_data$neg_log10_FDR)

sum(!is.finite(stage_plot_data$neg_log10_FDR))
# Top protein-coding by omnibus evidence
head(
  stage_plot_data[
    stage_plot_data$gene_type == "protein_coding",
  ][
    order(
      stage_plot_data[
        stage_plot_data$gene_type == "protein_coding",
      ]$adj.P.Val
    ),
    c(
      "gene_name",
      "max_abs_stage_effect",
      "AveExpr",
      "adj.P.Val"
    )
  ],
  15
)
# Top protein-coding by effect magnitude,
# restricted to omnibus-significant genes

tmp <- stage_plot_data[
  stage_plot_data$gene_type == "protein_coding" &
    stage_plot_data$adj.P.Val < 0.05,
]

head(
  tmp[
    order(tmp$max_abs_stage_effect, decreasing = TRUE),
    c(
      "gene_name",
      "max_abs_stage_effect",
      "AveExpr",
      "adj.P.Val"
    )
  ],
  15
)
# ============================================================
# S7.2 — Final Stage omnibus association overview
# ============================================================

library(ggplot2)
library(ggrepel)

stage_label_genes <- c(
  "ACKR1",
  "C7",
  "MMRN1",
  "ADH1B",
  "PIR",
  "MAGEA3",
  "ADIPOQ",
  "RCAN3"
)

stage_label_data <- stage_plot_data[
  stage_plot_data$gene_name %in% stage_label_genes,
]

stage_label_data <- stage_label_data[
  match(stage_label_genes, stage_label_data$gene_name),
]

# Separate RCAN3 so its label can be positioned independently
stage_label_main <- subset(
  stage_label_data,
  gene_name != "RCAN3"
)

stage_label_rcan3 <- subset(
  stage_label_data,
  gene_name == "RCAN3"
)

stage_overview_plot <- ggplot(
  stage_plot_data,
  aes(
    x = max_abs_stage_effect,
    y = neg_log10_FDR
  )
) +
  
  geom_point(
    data = subset(
      stage_plot_data,
      association_class == "Not significant"
    ),
    color = "grey78",
    alpha = 0.45,
    size = 1.25
  ) +
  
  geom_point(
    data = subset(
      stage_plot_data,
      association_class == "Stage-associated"
    ),
    color = "#7B4AB5",
    alpha = 0.75,
    size = 1.5
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +
  
  # Main labels
  geom_text_repel(
    data = stage_label_main,
    aes(label = gene_name),
    size = 3.4,
    seed = 42,
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  
  # RCAN3 positioned separately
  geom_text_repel(
    data = stage_label_rcan3,
    aes(label = gene_name),
    size = 3.4,
    seed = 42,
    nudge_x = -0.20,
    nudge_y = 0.12,
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  
  labs(
    title = "Transcriptome-wide association with pathological stage",
    subtitle = "TCGA-BRCA primary tumors (n = 1,063; 31,157 genes tested)",
    x = "Maximum absolute Stage-vs-I coefficient (voom log2-CPM)",
    y = expression(-log[10]("BH-FDR")),
    caption = "Dashed line: BH-FDR = 0.05. Effect magnitude is unsigned."
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"),
    plot.caption = element_text(
      size = 9,
      color = "grey35",
      hjust = 0
    ),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(10, 20, 10, 10)
  )

stage_overview_plot
# Export the final RCAN3-refined plot, not its superseded preview.
ggsave(
  file.path(output_dir, "figures/stage_omnibus_association_overview.png"),
  stage_overview_plot,
  width = 9,
  height = 6.8,
  dpi = 300
)

ggsave(
  file.path(output_dir, "figures/stage_omnibus_association_overview.pdf"),
  stage_overview_plot,
  width = 9,
  height = 6.8
)
# ============================================================
# S7.3 — Prepare patient-level Stage-expression data
# ============================================================

stage_selected_genes <- c(
  "ACKR1",
  "ADH1B",
  "PIR",
  "ADIPOQ"
)

stage_selected_info <- stage_results_annotated[
  stage_results_annotated$gene_name %in% stage_selected_genes,
]

stage_selected_info <- stage_selected_info[
  match(
    stage_selected_genes,
    stage_selected_info$gene_name
  ),
]

stage_selected_info[
  ,
  c(
    "ensembl_gene_id",
    "gene_name",
    "stage_groupStage.II",
    "stage_groupStage.III",
    "stage_groupStage.IV",
    "AveExpr",
    "F",
    "adj.P.Val"
  )
]
stage_selected_expression <- voom_stage_final$E[
  stage_selected_info$ensembl_gene_id,
  ,
  drop = FALSE
]

dim(stage_selected_expression)

identical(
  colnames(stage_selected_expression),
  stage_meta_final$patient_id
)
stage_expression_long <- data.frame(
  patient_id = rep(
    colnames(stage_selected_expression),
    times = nrow(stage_selected_expression)
  ),
  
  gene_name = rep(
    stage_selected_info$gene_name,
    each = ncol(stage_selected_expression)
  ),
  
  expression = as.vector(
    t(stage_selected_expression)
  )
)
stage_expression_long$stage_group <-
  stage_meta_final$stage_group[
    match(
      stage_expression_long$patient_id,
      stage_meta_final$patient_id
    )
  ]

stage_expression_long$stage_group <- factor(
  stage_expression_long$stage_group,
  levels = c(
    "Stage I",
    "Stage II",
    "Stage III",
    "Stage IV"
  )
)
dim(stage_expression_long)

sum(is.na(stage_expression_long$expression))
sum(is.na(stage_expression_long$stage_group))

table(stage_expression_long$stage_group)

table(
  stage_expression_long$gene_name,
  stage_expression_long$stage_group
)
# ============================================================
# S7.4 — Statistics labels for selected genes
# ============================================================

stage_selected_labels <- data.frame(
  gene_name = stage_selected_info$gene_name,
  label = paste0(
    "Omnibus FDR = ",
    formatC(
      stage_selected_info$adj.P.Val,
      format = "e",
      digits = 2
    )
  )
)
#rebuild
stage_expression_plot <- ggplot(
  stage_expression_long,
  aes(
    x = stage_group,
    y = expression
  )
) +
  
  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    linewidth = 0.45,
    fill = "grey92"
  ) +
  
  geom_jitter(
    width = 0.16,
    height = 0,
    alpha = 0.20,
    size = 0.85
  ) +
  
  geom_text(
    data = stage_selected_labels,
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.35,
    size = 3.1
  ) +
  
  facet_wrap(
    ~ gene_name,
    ncol = 2,
    scales = "free_y"
  ) +
  
  labs(
    title = "Expression of selected Stage-associated genes",
    subtitle = "TCGA-BRCA primary tumors (n = 1,063)",
    x = "Pathological stage",
    y = "Voom log2-CPM expression",
    caption = paste0(
      "Stage I n = 179; Stage II n = 615; ",
      "Stage III n = 249; Stage IV n = 20."
    )
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"),
    
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    
    axis.title = element_text(face = "bold"),
    
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    
    plot.caption = element_text(
      size = 9,
      color = "grey35",
      hjust = 0
    ),
    
    panel.spacing = grid::unit(1.1, "lines")
  )

stage_expression_plot
ggsave(
  file.path(output_dir, "figures/stage_selected_gene_expression.png"),
  stage_expression_plot,
  width = 9,
  height = 7.5,
  dpi = 300
)

ggsave(
  file.path(output_dir, "figures/stage_selected_gene_expression.pdf"),
  stage_expression_plot,
  width = 9,
  height = 7.5
)
# ============================================================
# S7.5 — Prepare Stage coefficient-pattern heatmap
# ============================================================

stage_pc_sig <- stage_results_annotated[
  stage_results_annotated$gene_type == "protein_coding" &
    stage_results_annotated$adj.P.Val < 0.05,
]

stage_pc_sig <- stage_pc_sig[
  order(stage_pc_sig$adj.P.Val),
]

stage_heatmap_top30 <- head(
  stage_pc_sig,
  30
)

stage_heatmap_matrix <- as.matrix(
  stage_heatmap_top30[
    ,
    c(
      "stage_groupStage.II",
      "stage_groupStage.III",
      "stage_groupStage.IV"
    )
  ]
)

rownames(stage_heatmap_matrix) <-
  stage_heatmap_top30$gene_name

colnames(stage_heatmap_matrix) <- c(
  "Stage II - I",
  "Stage III - I",
  "Stage IV - I"
)
dim(stage_heatmap_matrix)

sum(is.na(stage_heatmap_matrix))
sum(!is.finite(stage_heatmap_matrix))

rownames(stage_heatmap_matrix)

range(stage_heatmap_matrix)

stage_heatmap_matrix
# ============================================================
# S7.6 — Prepare long-format coefficient data
# ============================================================

stage_heatmap_long <- data.frame(
  gene_name = rep(
    rownames(stage_heatmap_matrix),
    times = ncol(stage_heatmap_matrix)
  ),
  
  comparison = rep(
    colnames(stage_heatmap_matrix),
    each = nrow(stage_heatmap_matrix)
  ),
  
  coefficient = as.vector(stage_heatmap_matrix)
)

stage_heatmap_long$gene_name <- factor(
  stage_heatmap_long$gene_name,
  levels = rev(rownames(stage_heatmap_matrix))
)

stage_heatmap_long$comparison <- factor(
  stage_heatmap_long$comparison,
  levels = c(
    "Stage II - I",
    "Stage III - I",
    "Stage IV - I"
  )
)

dim(stage_heatmap_long)

range(stage_heatmap_long$coefficient)
# ============================================================
# Figure 6 — Stage coefficient-pattern heatmap
# ============================================================

library(scales)

stage_coefficient_heatmap <- ggplot(
  stage_heatmap_long,
  aes(
    x = comparison,
    y = gene_name,
    fill = coefficient
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.35
  ) +
  
  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "white",
    high = "#D6534D",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    name = "Coefficient\n(log2-CPM)"
  ) +
  
  labs(
    title = "Stage-specific expression patterns among top Stage-associated genes",
    subtitle = "Top 30 protein-coding genes ranked by omnibus BH-FDR",
    x = NULL,
    y = NULL,
    caption = paste0(
      "Coefficients represent expression differences relative to Stage I. ",
      "Color scale capped at ±2.5 for visualization."
    )
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    
    plot.subtitle = element_text(
      color = "grey30"
    ),
    
    axis.text.x = element_text(
      face = "bold",
      angle = 0
    ),
    
    axis.text.y = element_text(
      size = 9
    ),
    
    panel.grid = element_blank(),
    
    legend.title = element_text(
      face = "bold"
    ),
    
    plot.caption = element_text(
      size = 9,
      color = "grey35",
      hjust = 0
    )
  )

stage_coefficient_heatmap
# CANONICAL final heatmap layout; one PNG/PDF export.
stage_coefficient_heatmap_final <-
  stage_coefficient_heatmap +
  
  labs(
    title = "Stage-specific expression patterns",
    
    subtitle = paste0(
      "Top 30 Stage-associated protein-coding genes\n",
      "ranked by omnibus BH-FDR"
    ),
    
    x = NULL,
    y = NULL,
    
    caption = paste0(
      "Coefficients represent expression differences relative to Stage I.\n",
      "Color scale capped at ±2.5 for visualization."
    )
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    
    plot.subtitle = element_text(
      size = 12,
      color = "grey30",
      lineheight = 1.05
    ),
    
    plot.caption = element_text(
      size = 9,
      color = "grey35",
      hjust = 0,
      lineheight = 1.05
    ),
    
    plot.margin = margin(
      t = 10,
      r = 20,
      b = 10,
      l = 10
    )
  )

stage_coefficient_heatmap_final

ggsave(
  file.path(output_dir, "figures/stage_coefficient_pattern_heatmap.png"),
  stage_coefficient_heatmap_final,
  width = 10,
  height = 8.5,
  dpi = 300
)

ggsave(
  file.path(output_dir, "figures/stage_coefficient_pattern_heatmap.pdf"),
  stage_coefficient_heatmap_final,
  width = 10,
  height = 8.5
)
# ============================================================
# S9.1 — Define common tested gene universe
# ============================================================

dim(age_results_annotated)
dim(stage_results_annotated)

names(age_results_annotated)
names(stage_results_annotated)
age_tested_genes <- age_results_annotated$ensembl_gene_id
stage_tested_genes <- stage_results_annotated$ensembl_gene_id

length(age_tested_genes)
length(stage_tested_genes)

length(unique(age_tested_genes))
length(unique(stage_tested_genes))

sum(is.na(age_tested_genes))
sum(is.na(stage_tested_genes))

anyDuplicated(age_tested_genes)
anyDuplicated(stage_tested_genes)
common_tested_genes <- intersect(
  age_tested_genes,
  stage_tested_genes
)

length(common_tested_genes)
all(
  age_tested_genes %in% stage_tested_genes
)

setdiff(
  age_tested_genes,
  stage_tested_genes
) |> length()

setdiff(
  stage_tested_genes,
  age_tested_genes
) |> length()
gene_universe_summary <- data.frame(
  universe = c(
    "Age tested",
    "Stage tested",
    "Tested in both",
    "Age only",
    "Stage only"
  ),
  n_genes = c(
    length(age_tested_genes),
    length(stage_tested_genes),
    length(common_tested_genes),
    length(setdiff(age_tested_genes, stage_tested_genes)),
    length(setdiff(stage_tested_genes, age_tested_genes))
  )
)

gene_universe_summary
age_integrated <- age_results_annotated[
  match(
    common_tested_genes,
    age_results_annotated$ensembl_gene_id
  ),
  c(
    "ensembl_gene_id",
    "gene_name",
    "gene_type",
    "beta_age_per_10y",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "FDR_significant"
  )
]

names(age_integrated) <- c(
  "ensembl_gene_id",
  "gene_name",
  "gene_type",
  "age_beta_10y",
  "age_AveExpr",
  "age_t",
  "age_P",
  "age_FDR",
  "age_significant"
)
stage_integrated <- stage_results_annotated[
  match(
    common_tested_genes,
    stage_results_annotated$ensembl_gene_id
  ),
  c(
    "ensembl_gene_id",
    "gene_name",
    "gene_type",
    "stage_groupStage.II",
    "stage_groupStage.III",
    "stage_groupStage.IV",
    "AveExpr",
    "F",
    "P.Value",
    "adj.P.Val",
    "FDR_significant"
  )
]

names(stage_integrated) <- c(
  "ensembl_gene_id",
  "stage_gene_name",
  "stage_gene_type",
  "stage_II_vs_I",
  "stage_III_vs_I",
  "stage_IV_vs_I",
  "stage_AveExpr",
  "stage_F",
  "stage_P",
  "stage_FDR",
  "stage_significant"
)
identical(
  age_integrated$ensembl_gene_id,
  stage_integrated$ensembl_gene_id
)

sum(
  age_integrated$gene_name != stage_integrated$stage_gene_name,
  na.rm = TRUE
)

sum(
  age_integrated$gene_type != stage_integrated$stage_gene_type,
  na.rm = TRUE
)

sum(is.na(age_integrated$gene_name))
sum(is.na(stage_integrated$stage_gene_name))
integrated_results <- data.frame(
  ensembl_gene_id = age_integrated$ensembl_gene_id,
  gene_name       = age_integrated$gene_name,
  gene_type       = age_integrated$gene_type,
  
  age_beta_10y    = age_integrated$age_beta_10y,
  age_AveExpr     = age_integrated$age_AveExpr,
  age_t           = age_integrated$age_t,
  age_P           = age_integrated$age_P,
  age_FDR         = age_integrated$age_FDR,
  age_significant = age_integrated$age_significant,
  
  stage_II_vs_I   = stage_integrated$stage_II_vs_I,
  stage_III_vs_I  = stage_integrated$stage_III_vs_I,
  stage_IV_vs_I   = stage_integrated$stage_IV_vs_I,
  stage_AveExpr   = stage_integrated$stage_AveExpr,
  stage_F         = stage_integrated$stage_F,
  stage_P         = stage_integrated$stage_P,
  stage_FDR       = stage_integrated$stage_FDR,
  stage_significant = stage_integrated$stage_significant,
  
  stringsAsFactors = FALSE
)
dim(integrated_results)

anyDuplicated(integrated_results$ensembl_gene_id)

sum(is.na(integrated_results$age_FDR))
sum(is.na(integrated_results$stage_FDR))

table(integrated_results$age_significant)
table(integrated_results$stage_significant)

summary(integrated_results$age_beta_10y)

summary(integrated_results$stage_FDR)
identical(
  age_integrated$ensembl_gene_id,
  stage_integrated$ensembl_gene_id
)

sum(
  age_integrated$gene_name != stage_integrated$stage_gene_name,
  na.rm = TRUE
)

sum(
  age_integrated$gene_type != stage_integrated$stage_gene_type,
  na.rm = TRUE
)

dim(integrated_results)

anyDuplicated(integrated_results$ensembl_gene_id)

table(integrated_results$age_significant)
table(integrated_results$stage_significant)
integrated_results$association_group <- with(
  integrated_results,
  ifelse(
    age_significant & stage_significant,
    "Both",
    ifelse(
      age_significant & !stage_significant,
      "Age only",
      ifelse(
        !age_significant & stage_significant,
        "Stage only",
        "Neither"
      )
    )
  )
)

integrated_results$association_group <- factor(
  integrated_results$association_group,
  levels = c(
    "Age only",
    "Stage only",
    "Both",
    "Neither"
  )
)

association_counts <- table(
  integrated_results$association_group
)

association_counts
sum(association_counts)

stopifnot(
  sum(association_counts) == 25551
)

association_summary <- data.frame(
  association_group = names(association_counts),
  n_genes = as.integer(association_counts),
  percent_common_universe =
    as.integer(association_counts) /
    nrow(integrated_results) * 100
)

association_summary
both_n <- sum(
  integrated_results$age_significant &
    integrated_results$stage_significant
)

age_sig_n <- sum(
  integrated_results$age_significant
)

stage_sig_n <- sum(
  integrated_results$stage_significant
)

both_n / age_sig_n * 100
both_n / stage_sig_n * 100
age_stage_table <- table(
  Age = integrated_results$age_significant,
  Stage = integrated_results$stage_significant
)

age_stage_table
addmargins(age_stage_table)
# ============================================================
# S9.3b — Expected overlap under independence
# ============================================================

expected_both <- (
  age_sig_n *
    stage_sig_n
) / nrow(integrated_results)

observed_both <- both_n

overlap_enrichment <- observed_both / expected_both

c(
  observed_both = observed_both,
  expected_both = expected_both,
  observed_expected_ratio = overlap_enrichment
)
age_stage_fisher <- fisher.test(
  age_stage_table
)

age_stage_fisher
age_stage_fisher$estimate
age_stage_fisher$conf.int
age_stage_fisher$p.value

# ============================================================
# S9.4a — Shared Age + Stage associated genes
# ============================================================

shared_results <- integrated_results[
  integrated_results$association_group == "Both",
]

dim(shared_results)

stopifnot(nrow(shared_results) == 513)
shared_results$age_direction <- ifelse(
  shared_results$age_beta_10y > 0,
  "Positive",
  "Negative"
)

table(shared_results$age_direction)

prop.table(
  table(shared_results$age_direction)
) * 100
shared_results$stage_pattern <- apply(
  shared_results[
    c(
      "stage_II_vs_I",
      "stage_III_vs_I",
      "stage_IV_vs_I"
    )
  ],
  1,
  function(x) {
    
    if (all(x > 0)) {
      "All positive"
    } else if (all(x < 0)) {
      "All negative"
    } else {
      "Mixed"
    }
    
  }
)

table(shared_results$stage_pattern)

prop.table(
  table(shared_results$stage_pattern)
) * 100
age_stage_pattern_table <- table(
  Age_direction = shared_results$age_direction,
  Stage_pattern = shared_results$stage_pattern
)

age_stage_pattern_table

addmargins(age_stage_pattern_table)
shared_gene_type_counts <- sort(
  table(shared_results$gene_type),
  decreasing = TRUE
)

shared_gene_type_counts

head(
  shared_gene_type_counts,
  15
)
head(
  prop.table(shared_gene_type_counts) * 100,
  15
)
shared_by_age <- shared_results[
  order(shared_results$age_FDR),
  c(
    "ensembl_gene_id",
    "gene_name",
    "gene_type",
    "age_beta_10y",
    "age_FDR",
    "stage_II_vs_I",
    "stage_III_vs_I",
    "stage_IV_vs_I",
    "stage_FDR"
  )
]

head(shared_by_age, 15)
shared_by_stage <- shared_results[
  order(shared_results$stage_FDR),
  c(
    "ensembl_gene_id",
    "gene_name",
    "gene_type",
    "age_beta_10y",
    "age_FDR",
    "stage_II_vs_I",
    "stage_III_vs_I",
    "stage_IV_vs_I",
    "stage_FDR"
  )
]

head(shared_by_stage, 15)
# ============================================================
# S9.4b — Gene-biotype representation in the shared set
# ============================================================

common_type_counts <- table(
  integrated_results$gene_type
)

shared_type_counts <- table(
  shared_results$gene_type
)

biotype_shared_summary <- data.frame(
  gene_type = names(common_type_counts),
  n_common = as.integer(common_type_counts),
  stringsAsFactors = FALSE
)

biotype_shared_summary$n_shared <- as.integer(
  shared_type_counts[
    match(
      biotype_shared_summary$gene_type,
      names(shared_type_counts)
    )
  ]
)

biotype_shared_summary$n_shared[
  is.na(biotype_shared_summary$n_shared)
] <- 0

biotype_shared_summary$percent_shared_within_type <-
  biotype_shared_summary$n_shared /
  biotype_shared_summary$n_common * 100

biotype_shared_summary <- biotype_shared_summary[
  order(
    biotype_shared_summary$n_shared,
    decreasing = TRUE
  ),
]

head(biotype_shared_summary, 15)
shared_results$stage_max_abs_effect <- apply(
  abs(
    shared_results[
      c(
        "stage_II_vs_I",
        "stage_III_vs_I",
        "stage_IV_vs_I"
      )
    ]
  ),
  1,
  max
)

summary(shared_results$age_beta_10y)

summary(abs(shared_results$age_beta_10y))

summary(shared_results$stage_max_abs_effect)
shared_top_effect <- shared_results[
  order(
    shared_results$stage_max_abs_effect,
    decreasing = TRUE
  ),
  c(
    "gene_name",
    "gene_type",
    "age_beta_10y",
    "age_FDR",
    "stage_II_vs_I",
    "stage_III_vs_I",
    "stage_IV_vs_I",
    "stage_max_abs_effect",
    "stage_FDR",
    "age_direction",
    "stage_pattern"
  )
]

head(shared_top_effect, 15)

# ============================================================
pipeline_event(ctx, "joint_inference")
# S9.5a — Joint cohort using canonical patient IDs
# ============================================================

head(age_meta_final$patient_id)
head(stage_meta_final$patient_id)

head(colnames(raw_count_matrix_final))

anyNA(age_meta_final$patient_id)
anyNA(stage_meta_final$patient_id)

anyDuplicated(age_meta_final$patient_id)
anyDuplicated(stage_meta_final$patient_id)

length(unique(age_meta_final$patient_id))
length(unique(stage_meta_final$patient_id))
identical(
  age_meta_final$patient_id,
  colnames(voom_age_final$E)
)

identical(
  stage_meta_final$patient_id,
  colnames(voom_stage_final$E)
)
all(
  age_meta_final$patient_id %in%
    colnames(raw_count_matrix_final)
)

all(
  stage_meta_final$patient_id %in%
    colnames(raw_count_matrix_final)
)
joint_patient_ids <- intersect(
  age_meta_final$patient_id,
  stage_meta_final$patient_id
)

length(joint_patient_ids)

length(setdiff(
  age_meta_final$patient_id,
  stage_meta_final$patient_id
))

length(setdiff(
  stage_meta_final$patient_id,
  age_meta_final$patient_id
))
joint_age_meta <- age_meta_final[
  match(
    joint_patient_ids,
    age_meta_final$patient_id
  ),
  ,
  drop = FALSE
]

joint_stage_meta <- stage_meta_final[
  match(
    joint_patient_ids,
    stage_meta_final$patient_id
  ),
  ,
  drop = FALSE
]

identical(
  joint_age_meta$patient_id,
  joint_stage_meta$patient_id
)

identical(
  joint_age_meta$patient_id,
  joint_patient_ids
)

anyDuplicated(joint_patient_ids)
all(
  joint_patient_ids %in%
    colnames(raw_count_matrix_final)
)

length(intersect(
  joint_patient_ids,
  colnames(raw_count_matrix_final)
))
identical(
  joint_age_meta$age_years,
  joint_stage_meta$age_years
)

identical(
  as.character(joint_age_meta$stage_group),
  as.character(joint_stage_meta$stage_group)
)
stopifnot(
  length(joint_patient_ids) > 5L, !anyNA(joint_patient_ids),
  !anyDuplicated(joint_patient_ids),
  identical(joint_age_meta$patient_id, joint_patient_ids),
  identical(joint_stage_meta$patient_id, joint_patient_ids),
  identical(joint_age_meta$age_years, joint_stage_meta$age_years),
  identical(as.character(joint_age_meta$stage_group),
            as.character(joint_stage_meta$stage_group)),
  all(joint_patient_ids %in% colnames(raw_count_matrix_final)),
  !any(final_qc_exclusions %in% joint_patient_ids),
  all(is.finite(joint_age_meta$age_years)), !anyNA(joint_stage_meta$stage_group)
)
# ============================================================
# S9.5b — Joint metadata and design matrix
# ============================================================

joint_meta <- joint_age_meta

# Explicit canonical row identifiers
rownames(joint_meta) <- joint_meta$patient_id

dim(joint_meta)
head(joint_meta[, c("patient_id", "age_years", "stage_group")])
sum(is.na(joint_meta$age_years))
sum(is.na(joint_meta$stage_group))

table(joint_meta$stage_group, useNA = "ifany")
summary(joint_meta$age_years)
joint_meta$age_10y <- joint_meta$age_years / 10
summary(joint_meta$age_10y)

all.equal(
  joint_meta$age_10y * 10,
  joint_meta$age_years
)
joint_meta$stage_group <- factor(
  as.character(joint_meta$stage_group),
  levels = c("Stage I", "Stage II", "Stage III", "Stage IV")
)

levels(joint_meta$stage_group)
table(joint_meta$stage_group)
design_joint <- model.matrix(
  ~ age_10y + stage_group,
  data = joint_meta,
  contrasts.arg = list(stage_group = contr.treatment(
    c("Stage I", "Stage II", "Stage III", "Stage IV"), base = 1L
  ))
)
stopifnot(
  identical(rownames(design_joint), joint_meta$patient_id),
  nrow(design_joint) == nrow(joint_meta),
  qr(design_joint)$rank == ncol(design_joint)
)

dim(design_joint)
colnames(design_joint)

head(design_joint)
qr(design_joint)$rank
ncol(design_joint)

qr(design_joint)$rank == ncol(design_joint)
identical(
  rownames(design_joint),
  joint_meta$patient_id
)

all(
  joint_meta$patient_id %in%
    colnames(raw_count_matrix_final)
)
aggregate(
  age_years ~ stage_group,
  data = joint_meta,
  FUN = function(x) c(
    n = length(x),
    mean = mean(x),
    sd = sd(x),
    median = median(x),
    min = min(x),
    max = max(x)
  )
)
# ============================================================
# S9.5c — Joint-specific expression preprocessing
# ============================================================

raw_count_joint <- raw_count_matrix_final[
  ,
  joint_meta$patient_id,
  drop = FALSE
]

dim(raw_count_joint)

identical(
  colnames(raw_count_joint),
  joint_meta$patient_id
)

identical(
  colnames(raw_count_joint),
  rownames(design_joint)
)

anyNA(raw_count_joint)
anyDuplicated(colnames(raw_count_joint))
all(raw_count_joint >= 0)

all(
  raw_count_joint ==
    round(raw_count_joint)
)
library(edgeR)
library(limma)

dge_joint <- DGEList(
  counts = raw_count_joint
)

dim(dge_joint)

summary(dge_joint$samples$lib.size)
keep_joint <- filterByExpr(
  dge_joint,
  design = design_joint
)

table(keep_joint)

sum(keep_joint)
dge_joint_filtered <- dge_joint[
  keep_joint,
  ,
  keep.lib.sizes = FALSE
]

dim(dge_joint_filtered)
dge_joint_norm <- calcNormFactors(
  dge_joint_filtered,
  method = "TMM"
)

summary(dge_joint_norm$samples$norm.factors)

dge_joint_norm$samples
effective_lib_sizes_joint <-
  dge_joint_norm$samples$lib.size *
  dge_joint_norm$samples$norm.factors

summary(effective_lib_sizes_joint)
voom_joint <- voom(
  dge_joint_norm,
  design_joint,
  plot = TRUE
)
dim(voom_joint$E)
dim(voom_joint$weights)

identical(
  colnames(voom_joint$E),
  rownames(design_joint)
)

identical(
  rownames(voom_joint$E),
  rownames(dge_joint_norm$counts)
)

anyNA(voom_joint$E)
anyNA(voom_joint$weights)

range(voom_joint$weights)
stopifnot(
  identical(colnames(voom_joint$E), joint_meta$patient_id),
  identical(dim(voom_joint$E), dim(voom_joint$weights)),
  identical(rownames(voom_joint$E), rownames(dge_joint_norm$counts)),
  all(is.finite(voom_joint$E)),
  all(is.finite(voom_joint$weights)), all(voom_joint$weights > 0),
  all(is.finite(dge_joint_norm$samples$norm.factors)),
  all(dge_joint_norm$samples$norm.factors > 0)
)
# ============================================================
# S9.5d — Joint Age + Stage inference
# ============================================================

fit_joint <- lmFit(
  voom_joint,
  design_joint
)
pipeline_assert_fit(fit_joint, voom_joint, design_joint)

fit_joint_eb <- eBayes(fit_joint)

dim(fit_joint$coefficients)

colnames(fit_joint$coefficients)

summary(fit_joint$df.residual)
joint_age_results <- topTable(
  fit_joint_eb,
  coef = "age_10y",
  number = Inf,
  sort.by = "none"
)
pipeline_assert_results(joint_age_results)

dim(joint_age_results)

head(joint_age_results)

sum(joint_age_results$adj.P.Val < 0.05)

sum(joint_age_results$adj.P.Val < 0.01)

sum(joint_age_results$adj.P.Val < 0.001)

summary(joint_age_results$logFC)
table(
  Significant = joint_age_results$adj.P.Val < 0.05,
  Direction = ifelse(
    joint_age_results$logFC > 0,
    "Positive",
    "Negative"
  )
)
joint_stage_coef_names <- c(
  "stage_groupStage II",
  "stage_groupStage III",
  "stage_groupStage IV"
)

joint_stage_omnibus <- topTable(
  fit_joint_eb,
  coef = joint_stage_coef_names,
  number = Inf,
  sort.by = "none"
)
pipeline_assert_results(joint_stage_omnibus)

dim(joint_stage_omnibus)

head(joint_stage_omnibus)

sum(joint_stage_omnibus$adj.P.Val < 0.05)

sum(joint_stage_omnibus$adj.P.Val < 0.01)

sum(joint_stage_omnibus$adj.P.Val < 0.001)

summary(joint_stage_omnibus$F)
summary(joint_stage_omnibus$adj.P.Val)
identical(
  rownames(joint_age_results),
  rownames(voom_joint$E)
)

identical(
  rownames(joint_stage_omnibus),
  rownames(voom_joint$E)
)

identical(
  rownames(joint_age_results),
  rownames(joint_stage_omnibus)
)
# ============================================================
# S9.6a — Primary vs Stage-adjusted Age associations
# ============================================================

age_joint_common_genes <- intersect(
  rownames(age_results_final),
  rownames(joint_age_results)
)

length(age_joint_common_genes)

length(setdiff(
  rownames(age_results_final),
  rownames(joint_age_results)
))

length(setdiff(
  rownames(joint_age_results),
  rownames(age_results_final)
))
age_primary_matched <- age_results_final[
  match(
    age_joint_common_genes,
    rownames(age_results_final)
  ),
  ,
  drop = FALSE
]

age_adjusted_matched <- joint_age_results[
  match(
    age_joint_common_genes,
    rownames(joint_age_results)
  ),
  ,
  drop = FALSE
]

identical(
  rownames(age_primary_matched),
  rownames(age_adjusted_matched)
)
age_primary_beta_10y <-
  age_primary_matched$beta_age_per_10y
age_adjusted_beta_10y <-
  age_adjusted_matched$logFC
age_integrated_check <- integrated_results[
  match(age_joint_common_genes, integrated_results$ensembl_gene_id), , drop = FALSE
]
stopifnot(
  identical(age_integrated_check$ensembl_gene_id, age_joint_common_genes),
  length(age_primary_beta_10y) == length(age_adjusted_beta_10y),
  all(is.finite(age_primary_beta_10y)), all(is.finite(age_adjusted_beta_10y)),
  isTRUE(all.equal(age_primary_beta_10y, age_integrated_check$age_beta_10y,
                   tolerance = 1e-12))
)
age_beta_change <- age_adjusted_beta_10y - age_primary_beta_10y

age_primary_sig <-
  age_primary_matched$adj.P.Val < 0.05

age_adjusted_sig <-
  age_adjusted_matched$adj.P.Val < 0.05

age_sig_stability <- table(
  Primary = age_primary_sig,
  Adjusted = age_adjusted_sig
)

age_sig_stability

addmargins(age_sig_stability)
age_retention_rate <-
  sum(age_primary_sig & age_adjusted_sig) /
  sum(age_primary_sig) * 100

age_retention_rate

length(age_primary_beta_10y)
length(age_adjusted_beta_10y)

summary(age_primary_beta_10y)
summary(age_adjusted_beta_10y)
cor(
  age_primary_beta_10y,
  age_adjusted_beta_10y,
  method = "pearson"
)

cor(
  age_primary_beta_10y,
  age_adjusted_beta_10y,
  method = "spearman"
)

mean(
  sign(age_primary_beta_10y) ==
    sign(age_adjusted_beta_10y)
) * 100

summary(age_beta_change)

summary(abs(age_beta_change))

quantile(
  abs(age_beta_change),
  probs = c(0.90, 0.95, 0.99)
)
mean(
  abs(age_beta_change) < 0.05
) * 100

mean(
  abs(age_beta_change) < 0.10
) * 100
# ============================================================
# S9.6b — Primary vs Age-adjusted Stage associations
# ============================================================

stage_joint_common_genes <- intersect(
  rownames(stage_omnibus),
  rownames(joint_stage_omnibus)
)

length(stage_joint_common_genes)

length(setdiff(
  rownames(stage_omnibus),
  rownames(joint_stage_omnibus)
))

length(setdiff(
  rownames(joint_stage_omnibus),
  rownames(stage_omnibus)
))
stage_primary_matched <- stage_omnibus[
  match(
    stage_joint_common_genes,
    rownames(stage_omnibus)
  ),
  ,
  drop = FALSE
]

stage_adjusted_matched <- joint_stage_omnibus[
  match(
    stage_joint_common_genes,
    rownames(joint_stage_omnibus)
  ),
  ,
  drop = FALSE
]

identical(
  rownames(stage_primary_matched),
  rownames(stage_adjusted_matched)
)
stage_primary_sig <-
  stage_primary_matched$adj.P.Val < 0.05

stage_adjusted_sig <-
  stage_adjusted_matched$adj.P.Val < 0.05

stage_sig_stability <- table(
  Primary = stage_primary_sig,
  Adjusted = stage_adjusted_sig
)

stage_sig_stability
addmargins(stage_sig_stability)

stage_retention_rate <-
  sum(stage_primary_sig & stage_adjusted_sig) /
  sum(stage_primary_sig) * 100

stage_retention_rate
names(stage_omnibus)
names(joint_stage_omnibus)
stage_coef_names <- c(
  "stage_groupStage.II",
  "stage_groupStage.III",
  "stage_groupStage.IV"
)
stage_effect_stability <- do.call(
  rbind,
  lapply(stage_coef_names, function(coef_name) {
    
    primary_beta <- stage_primary_matched[[coef_name]]
    adjusted_beta <- stage_adjusted_matched[[coef_name]]
    
    beta_change <- adjusted_beta - primary_beta
    abs_change <- abs(beta_change)
    
    data.frame(
      contrast = coef_name,
      
      pearson = cor(
        primary_beta,
        adjusted_beta,
        method = "pearson"
      ),
      
      spearman = cor(
        primary_beta,
        adjusted_beta,
        method = "spearman"
      ),
      
      sign_concordance_pct =
        mean(
          sign(primary_beta) ==
            sign(adjusted_beta)
        ) * 100,
      
      median_abs_change =
        median(abs_change),
      
      mean_abs_change =
        mean(abs_change),
      
      p90_abs_change =
        unname(quantile(abs_change, 0.90)),
      
      p95_abs_change =
        unname(quantile(abs_change, 0.95)),
      
      p99_abs_change =
        unname(quantile(abs_change, 0.99)),
      
      max_abs_change =
        max(abs_change),
      
      pct_abs_change_lt_0.05 =
        mean(abs_change < 0.05) * 100,
      
      pct_abs_change_lt_0.10 =
        mean(abs_change < 0.10) * 100
    )
  })
)

stage_effect_stability
for (coef_name in stage_coef_names) {
  
  primary_beta <-
    stage_primary_matched[[coef_name]]
  
  adjusted_beta <-
    stage_adjusted_matched[[coef_name]]
  
  beta_change <-
    adjusted_beta - primary_beta
  
  cat("\n============================\n")
  cat(coef_name, "\n")
  cat("============================\n")
  
  cat("\nChange (adjusted - primary):\n")
  print(summary(beta_change))
  
  cat("\nAbsolute change:\n")
  print(summary(abs(beta_change)))
}
sapply(
  stage_coef_names,
  function(x) {
    sum(is.na(stage_primary_matched[[x]]))
  }
)

sapply(
  stage_coef_names,
  function(x) {
    sum(is.na(stage_adjusted_matched[[x]]))
  }
)
# ============================================================
# S9.7 — Adjusted status of the 513 primary shared genes
# ============================================================

shared_gene_ids <- shared_results$ensembl_gene_id

length(shared_gene_ids)
anyDuplicated(shared_gene_ids)

shared_age_adjusted <- joint_age_results[
  match(
    shared_gene_ids,
    rownames(joint_age_results)
  ),
  ,
  drop = FALSE
]

shared_stage_adjusted <- joint_stage_omnibus[
  match(
    shared_gene_ids,
    rownames(joint_stage_omnibus)
  ),
  ,
  drop = FALSE
]

sum(is.na(rownames(shared_age_adjusted)))
sum(is.na(rownames(shared_stage_adjusted)))

identical(
  rownames(shared_age_adjusted),
  shared_gene_ids
)

identical(
  rownames(shared_stage_adjusted),
  shared_gene_ids
)
shared_age_sig_adjusted <-
  shared_age_adjusted$adj.P.Val < 0.05

shared_stage_sig_adjusted <-
  shared_stage_adjusted$adj.P.Val < 0.05

shared_adjusted_status <- table(
  Age_adjusted = shared_age_sig_adjusted,
  Stage_adjusted = shared_stage_sig_adjusted
)

shared_adjusted_status
both_retained_n <- sum(
  shared_age_sig_adjusted &
    shared_stage_sig_adjusted
)

both_retained_pct <-
  both_retained_n /
  length(shared_gene_ids) * 100

c(
  shared_primary = length(shared_gene_ids),
  both_retained = both_retained_n,
  both_retained_pct = both_retained_pct
)

addmargins(shared_adjusted_status)
# ============================================================
pipeline_event(ctx, "integration")
# S9.8a — Step 2
# Build canonical final integrated master table
# ============================================================

master_gene_ids <- integrated_results$ensembl_gene_id

joint_age_master <- joint_age_results[
  match(
    master_gene_ids,
    rownames(joint_age_results)
  ),
  ,
  drop = FALSE
]

joint_stage_master <- joint_stage_omnibus[
  match(
    master_gene_ids,
    rownames(joint_stage_omnibus)
  ),
  ,
  drop = FALSE
]


# Alignment checks
identical(
  rownames(joint_age_master),
  master_gene_ids
)

identical(
  rownames(joint_stage_master),
  master_gene_ids
)

anyNA(joint_age_master)
anyNA(joint_stage_master)
integrated_master_final <- data.frame(
  
  # -----------------------------
  # Gene annotation
  # -----------------------------
  ensembl_gene_id =
    integrated_results$ensembl_gene_id,
  
  gene_name =
    integrated_results$gene_name,
  
  gene_type =
    integrated_results$gene_type,
  
  # -----------------------------
  # Primary Age analysis
  # -----------------------------
  age_primary_beta_10y =
    integrated_results$age_beta_10y,
  
  age_primary_AveExpr =
    integrated_results$age_AveExpr,
  
  age_primary_t =
    integrated_results$age_t,
  
  age_primary_P =
    integrated_results$age_P,
  
  age_primary_FDR =
    integrated_results$age_FDR,
  
  age_primary_significant =
    integrated_results$age_significant,
  
  # -----------------------------
  # Primary Stage analysis
  # -----------------------------
  stage_primary_II_vs_I =
    integrated_results$stage_II_vs_I,
  
  stage_primary_III_vs_I =
    integrated_results$stage_III_vs_I,
  
  stage_primary_IV_vs_I =
    integrated_results$stage_IV_vs_I,
  
  stage_primary_AveExpr =
    integrated_results$stage_AveExpr,
  
  stage_primary_F =
    integrated_results$stage_F,
  
  stage_primary_P =
    integrated_results$stage_P,
  
  stage_primary_FDR =
    integrated_results$stage_FDR,
  
  stage_primary_significant =
    integrated_results$stage_significant,
  
  # -----------------------------
  # Primary integrated category
  # -----------------------------
  association_group =
    integrated_results$association_group,
  
  # -----------------------------
  # Age adjusted for Stage
  # -----------------------------
  age_adjusted_beta_10y =
    joint_age_master$logFC,
  
  age_adjusted_AveExpr =
    joint_age_master$AveExpr,
  
  age_adjusted_t =
    joint_age_master$t,
  
  age_adjusted_P =
    joint_age_master$P.Value,
  
  age_adjusted_FDR =
    joint_age_master$adj.P.Val,
  
  age_adjusted_significant =
    joint_age_master$adj.P.Val < 0.05,
  
  # -----------------------------
  # Stage adjusted for Age
  # -----------------------------
  stage_adjusted_II_vs_I =
    joint_stage_master$stage_groupStage.II,
  
  stage_adjusted_III_vs_I =
    joint_stage_master$stage_groupStage.III,
  
  stage_adjusted_IV_vs_I =
    joint_stage_master$stage_groupStage.IV,
  
  stage_adjusted_AveExpr =
    joint_stage_master$AveExpr,
  
  stage_adjusted_F =
    joint_stage_master$F,
  
  stage_adjusted_P =
    joint_stage_master$P.Value,
  
  stage_adjusted_FDR =
    joint_stage_master$adj.P.Val,
  
  stage_adjusted_significant =
    joint_stage_master$adj.P.Val < 0.05,
  
  stringsAsFactors = FALSE
)
integrated_master_final$adjusted_association_group <- ifelse(
  integrated_master_final$age_adjusted_significant &
    integrated_master_final$stage_adjusted_significant,
  "Both",
  
  ifelse(
    integrated_master_final$age_adjusted_significant,
    "Age only",
    
    ifelse(
      integrated_master_final$stage_adjusted_significant,
      "Stage only",
      "Neither"
    )
  )
)
integrated_master_final$association_transition <- paste(
  integrated_master_final$association_group,
  "->",
  integrated_master_final$adjusted_association_group
)
pipeline_assert_ids(integrated_master_final$ensembl_gene_id, "integrated genes")
stopifnot(!anyNA(integrated_master_final[, setdiff(names(integrated_master_final), c("gene_name", "gene_type")), drop = FALSE]),
          all(integrated_master_final$age_primary_significant == (integrated_master_final$age_primary_FDR < 0.05)),
          all(integrated_master_final$stage_primary_significant == (integrated_master_final$stage_primary_FDR < 0.05)),
          all(integrated_master_final$age_adjusted_significant == (integrated_master_final$age_adjusted_FDR < 0.05)),
          all(integrated_master_final$stage_adjusted_significant == (integrated_master_final$stage_adjusted_FDR < 0.05)))
dim(integrated_master_final)

anyNA(integrated_master_final$ensembl_gene_id)
anyDuplicated(integrated_master_final$ensembl_gene_id)

identical(
  integrated_master_final$ensembl_gene_id,
  master_gene_ids
)

table(
  integrated_master_final$association_group,
  useNA = "ifany"
)

table(
  integrated_master_final$adjusted_association_group,
  useNA = "ifany"
)

table(
  Primary = integrated_master_final$association_group,
  Adjusted = integrated_master_final$adjusted_association_group
)

sum(integrated_master_final$age_primary_significant)
sum(integrated_master_final$stage_primary_significant)

sum(integrated_master_final$age_adjusted_significant)
sum(integrated_master_final$stage_adjusted_significant)
shared_master_check <- integrated_master_final[
  integrated_master_final$association_group == "Both",
]

nrow(shared_master_check)

table(
  shared_master_check$adjusted_association_group
)
# ============================================================
# S9.8b — Figure 7
# Step 1: Prepare primary association category counts
# ============================================================

figure7_group_order <- c(
  "Age only",
  "Stage only",
  "Both",
  "Neither"
)

figure7_counts <- as.data.frame(
  table(
    factor(
      integrated_master_final$association_group,
      levels = figure7_group_order
    )
  )
)

names(figure7_counts) <- c(
  "association_group",
  "n_genes"
)

figure7_counts$percentage <-
  100 * figure7_counts$n_genes /
  nrow(integrated_master_final)

figure7_counts$label <- sprintf(
  "%s (%.2f%%)",
  format(
    figure7_counts$n_genes,
    big.mark = ",",
    scientific = FALSE
  ),
  figure7_counts$percentage
)

figure7_counts
# Integrity checks

sum(figure7_counts$n_genes)

anyNA(figure7_counts)

identical(
  as.integer(figure7_counts$n_genes),
  c(8550L, 653L, 513L, 15835L)
)

sum(figure7_counts$percentage)
# ============================================================
# S9.8b — Figure 7
# Primary Age–Stage association categories
# ============================================================

library(ggplot2)

figure7_colors <- c(
  "Age only"   = "#4F6D9A",
  "Stage only" = "#D59A56",
  "Both"       = "#8664A8",
  "Neither"    = "#BFC5CD"
)

figure7_plot <- ggplot(
  figure7_counts,
  aes(
    x = n_genes,
    y = factor(
      association_group,
      levels = rev(figure7_group_order)
    ),
    fill = association_group
  )
) +
  geom_col(
    width = 0.64,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = label),
    hjust = -0.12,
    size = 3.8,
    color = "#30343B"
  ) +
  scale_fill_manual(
    values = figure7_colors
  ) +
  scale_x_continuous(
    labels = scales::label_comma(),
    limits = c(0, 19500),
    breaks = seq(0, 15000, by = 5000),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Age–Stage transcriptomic association categories",
    subtitle = "TCGA-BRCA | Common primary testing universe (n = 25,551 genes)",
    x = "Number of genes",
    y = NULL,
    caption = paste(
      "Primary Age and Stage analyses; BH-FDR < 0.05.",
      "Categories represent statistical associations, not causal specificity."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      color = "#252A33"
    ),
    plot.subtitle = element_text(
      size = 10,
      color = "#606873",
      margin = margin(b = 15)
    ),
    axis.text.y = element_text(
      size = 11,
      face = "bold",
      color = "#30343B"
    ),
    axis.text.x = element_text(
      size = 10,
      color = "#606873"
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(
      size = 8.5,
      color = "#69717D",
      hjust = 0,
      margin = margin(t = 15)
    ),
    plot.margin = margin(15, 25, 15, 15)
  )

print(figure7_plot)
# ============================================================
# S9.8c — Figure 8
# Step 1: Prepare matched primary vs adjusted coefficients
# ============================================================

figure8_data <- rbind(
  data.frame(
    gene_id = integrated_master_final$ensembl_gene_id,
    comparison = "Age (per 10 years)",
    primary = integrated_master_final$age_primary_beta_10y,
    adjusted = integrated_master_final$age_adjusted_beta_10y
  ),
  data.frame(
    gene_id = integrated_master_final$ensembl_gene_id,
    comparison = "Stage II vs I",
    primary = integrated_master_final$stage_primary_II_vs_I,
    adjusted = integrated_master_final$stage_adjusted_II_vs_I
  ),
  data.frame(
    gene_id = integrated_master_final$ensembl_gene_id,
    comparison = "Stage III vs I",
    primary = integrated_master_final$stage_primary_III_vs_I,
    adjusted = integrated_master_final$stage_adjusted_III_vs_I
  ),
  data.frame(
    gene_id = integrated_master_final$ensembl_gene_id,
    comparison = "Stage IV vs I",
    primary = integrated_master_final$stage_primary_IV_vs_I,
    adjusted = integrated_master_final$stage_adjusted_IV_vs_I
  )
)

figure8_data$comparison <- factor(
  figure8_data$comparison,
  levels = c(
    "Age (per 10 years)",
    "Stage II vs I",
    "Stage III vs I",
    "Stage IV vs I"
  )
)
dim(figure8_data)

table(figure8_data$comparison)

anyNA(figure8_data)

all(is.finite(figure8_data$primary))
all(is.finite(figure8_data$adjusted))

aggregate(
  cbind(primary, adjusted) ~ comparison,
  data = figure8_data,
  FUN = length
)
# ============================================================
# S9.8c — Figure 8
# Step 2: Quantify primary vs adjusted effect stability
# ============================================================

figure8_split <- split(
  figure8_data,
  figure8_data$comparison
)

figure8_metrics <- do.call(
  rbind,
  lapply(names(figure8_split), function(comp) {
    
    df <- figure8_split[[comp]]
    
    delta <- df$adjusted - df$primary
    
    data.frame(
      comparison = comp,
      n_genes = nrow(df),
      
      pearson_r = cor(
        df$primary,
        df$adjusted,
        method = "pearson"
      ),
      
      spearman_rho = cor(
        df$primary,
        df$adjusted,
        method = "spearman"
      ),
      
      sign_concordance_pct = 100 * mean(
        sign(df$primary) == sign(df$adjusted)
      ),
      
      median_abs_delta = median(abs(delta)),
      
      p95_abs_delta = as.numeric(
        quantile(abs(delta), 0.95)
      ),
      
      max_abs_delta = max(abs(delta))
    )
  })
)

rownames(figure8_metrics) <- NULL

figure8_metrics$comparison <- factor(
  figure8_metrics$comparison,
  levels = levels(figure8_data$comparison)
)

figure8_metrics <- figure8_metrics[
  order(figure8_metrics$comparison),
]

print(
  figure8_metrics,
  digits = 6,
  row.names = FALSE
)
# Exactly four comparisons
nrow(figure8_metrics)

# Each comparison contains the same number of genes
all(figure8_metrics$n_genes == 25551)

# Correlations within valid bounds
all(abs(figure8_metrics$pearson_r) <= 1)
all(abs(figure8_metrics$spearman_rho) <= 1)

# Non-negative absolute coefficient changes
all(figure8_metrics$median_abs_delta >= 0)

# Age Pearson correlation reproduces S9.6
all.equal(
  figure8_metrics$pearson_r[
    figure8_metrics$comparison == "Age (per 10 years)"
  ],
  0.997873,
  tolerance = 1e-5
)
# ============================================================
# S9.8c — Figure 8
# Step 3: Primary vs adjusted effect-estimate stability
# ============================================================

library(ggplot2)

# Panel annotations
figure8_metrics$annotation <- sprintf(
  "Pearson r = %.4f\nMedian |Delta beta| = %.4f",
  figure8_metrics$pearson_r,
  figure8_metrics$median_abs_delta
)

# Construct plot
figure8_plot <- ggplot(
  figure8_data,
  aes(
    x = primary,
    y = adjusted
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "#A0A5AE"
  ) +
  geom_point(
    color = "#526F9E",
    alpha = 0.16,
    size = 0.45
  ) +
  facet_wrap(
    ~ comparison,
    ncol = 2,
    scales = "free"
  ) +
  coord_equal() +
  labs(
    title = "Stability of Age and Stage effect estimates",
    subtitle = paste(
      "Primary versus mutually adjusted models |",
      "25,551 matched genes per comparison"
    ),
    x = "Primary model coefficient",
    y = "Mutually adjusted model coefficient",
    caption = paste0(
      "Dashed line indicates equality (y = x). ",
      "Age coefficients are per 10 years; Stage coefficients are relative to Stage I.\n",
      "Comparisons also reflect changes in cohort and preprocessing. ",
      "High concordance does not establish causal independence."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      color = "#252A33"
    ),
    plot.subtitle = element_text(
      color = "#606873",
      margin = margin(b = 12)
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(
      hjust = 0,
      size = 8,
      lineheight = 1.3,
      margin = margin(t = 12)
    ),
    plot.margin = margin(15, 20, 15, 15)
  )

# Defer rendering until the valid final coordinate system is installed below.
# The initial coord_equal/free-facet combination is incompatible with ggplot2 4.
# ============================================================
# S9.8c — Figure 8
# Fix: panel-specific symmetric axis limits
# ============================================================

# Determine a shared x/y range within each comparison
figure8_axis_limits <- do.call(
  rbind,
  lapply(split(figure8_data, figure8_data$comparison), function(df) {
    
    lim <- max(abs(c(df$primary, df$adjusted)))
    
    data.frame(
      comparison = as.character(df$comparison[1]),
      primary = c(-lim, lim),
      adjusted = c(-lim, lim)
    )
  })
)

figure8_axis_limits$comparison <- factor(
  figure8_axis_limits$comparison,
  levels = levels(figure8_data$comparison)
)

# Rebuild the plot with panel-specific axis ranges
figure8_plot <- figure8_plot +
  geom_blank(
    data = figure8_axis_limits,
    mapping = aes(x = primary, y = adjusted),
    inherit.aes = FALSE
  ) +
  coord_cartesian()

print(figure8_plot)
# ============================================================
# S9.8c — Figure 8
# Step 4: Statistical annotations and layout refinement
# ============================================================

figure8_annotation <- data.frame(
  comparison = figure8_metrics$comparison,
  label = sprintf(
    "Pearson r = %.4f\nMedian |Delta beta| = %.4f",
    figure8_metrics$pearson_r,
    figure8_metrics$median_abs_delta
  )
)

figure8_plot <- figure8_plot +
  geom_text(
    data = figure8_annotation,
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.08,
    vjust = 1.2,
    size = 3.1,
    lineheight = 1.3,
    color = "#30343B"
  ) +
  labs(
    subtitle = "Primary vs mutually adjusted models | 25,551 genes per panel",
    caption = paste0(
      "Dashed line: equality (y = x). ",
      "Age coefficients: per 10 years; Stage: relative to Stage I.\n",
      "Models also differ in cohort and preprocessing; ",
      "concordance does not establish causal independence."
    )
  ) +
  theme(
    plot.subtitle = element_text(
      size = 10,
      color = "#606873",
      margin = margin(b = 14)
    ),
    plot.caption = element_text(
      size = 8,
      hjust = 0,
      lineheight = 1.35,
      margin = margin(t = 16)
    ),
    panel.spacing = grid::unit(1.5, "lines"),
    plot.margin = margin(18, 25, 20, 18)
  )

print(figure8_plot)
# ============================================================
# S9.8d — Figure 9
# Step 1: Define stable shared-gene candidate pool
# ============================================================

figure9_candidates <- integrated_master_final[
  integrated_master_final$association_group == "Both" &
    integrated_master_final$adjusted_association_group == "Both",
  ,
  drop = FALSE
]

# Integrity checks
dim(figure9_candidates)

anyNA(figure9_candidates$ensembl_gene_id)

anyDuplicated(figure9_candidates$ensembl_gene_id)

all(
  figure9_candidates$age_primary_significant &
    figure9_candidates$stage_primary_significant &
    figure9_candidates$age_adjusted_significant &
    figure9_candidates$stage_adjusted_significant
)

summary(
  figure9_candidates[, c(
    "age_adjusted_FDR",
    "stage_adjusted_FDR"
  )]
)
# ============================================================
# S9.8d — Figure 9
# Step 2: Reproducible ranking of stable shared genes
# ============================================================

figure9_candidates$rank_age <- rank(
  figure9_candidates$age_adjusted_FDR,
  ties.method = "min"
)

figure9_candidates$rank_stage <- rank(
  figure9_candidates$stage_adjusted_FDR,
  ties.method = "min"
)

figure9_candidates$worst_rank <- pmax(
  figure9_candidates$rank_age,
  figure9_candidates$rank_stage
)

figure9_candidates$rank_sum <-
  figure9_candidates$rank_age +
  figure9_candidates$rank_stage

figure9_ranked <- figure9_candidates[
  order(
    figure9_candidates$worst_rank,
    figure9_candidates$rank_sum,
    figure9_candidates$ensembl_gene_id
  ),
  ,
  drop = FALSE
]

rownames(figure9_ranked) <- NULL

figure9_ranked$selection_rank <- seq_len(
  nrow(figure9_ranked)
)

figure9_top30 <- head(figure9_ranked, 30)
dim(figure9_ranked)
dim(figure9_top30)

anyDuplicated(figure9_top30$ensembl_gene_id)

all(
  figure9_top30$age_adjusted_FDR < 0.05 &
    figure9_top30$stage_adjusted_FDR < 0.05
)

all(diff(figure9_ranked$worst_rank) >= 0)
figure9_top30[
  ,
  c(
    "selection_rank",
    "ensembl_gene_id",
    "gene_name",
    "gene_type",
    "age_adjusted_FDR",
    "stage_adjusted_FDR",
    "rank_age",
    "rank_stage",
    "worst_rank"
  )
]
sum(is.na(figure9_top30$gene_name))

sum(
  is.na(figure9_top30$gene_name) |
    trimws(figure9_top30$gene_name) == ""
)
# ============================================================
# S9.8e — Figure 9
# Step 1: Build adjusted coefficient matrix
# ============================================================

figure9_coef_matrix <- as.matrix(
  figure9_top30[
    ,
    c(
      "age_adjusted_beta_10y",
      "stage_adjusted_II_vs_I",
      "stage_adjusted_III_vs_I",
      "stage_adjusted_IV_vs_I"
    )
  ]
)

colnames(figure9_coef_matrix) <- c(
  "Age / 10 years",
  "Stage II vs I",
  "Stage III vs I",
  "Stage IV vs I"
)

rownames(figure9_coef_matrix) <-
  figure9_top30$gene_name

# Integrity checks
dim(figure9_coef_matrix)

anyNA(figure9_coef_matrix)

all(is.finite(figure9_coef_matrix))

anyDuplicated(rownames(figure9_coef_matrix))

summary(as.vector(figure9_coef_matrix))

apply(
  figure9_coef_matrix,
  2,
  range
)
# ============================================================
# S9.8e — Figure 9
# Step 2: Zero-preserving column-wise scaling
# ============================================================

# Maximum absolute coefficient in each column
figure9_scale_factors <- apply(
  figure9_coef_matrix,
  2,
  function(x) max(abs(x))
)

figure9_scale_factors

# Ensure all scale factors are positive and finite
stopifnot(
  all(is.finite(figure9_scale_factors)),
  all(figure9_scale_factors > 0)
)

# Scale each column independently
figure9_scaled_matrix <- sweep(
  figure9_coef_matrix,
  2,
  figure9_scale_factors,
  FUN = "/"
)

# Validation
dim(figure9_scaled_matrix)

anyNA(figure9_scaled_matrix)

all(is.finite(figure9_scaled_matrix))

all(abs(figure9_scaled_matrix) <= 1 + 1e-12)

identical(
  sign(figure9_scaled_matrix),
  sign(figure9_coef_matrix)
)

apply(
  figure9_scaled_matrix,
  2,
  range
)

# Verify that scaling can be reversed
all.equal(
  sweep(
    figure9_scaled_matrix,
    2,
    figure9_scale_factors,
    FUN = "*"
  ),
  figure9_coef_matrix
)
# ============================================================
# S9.8e — Figure 9
# Step 3: Integrated adjusted-coefficient heatmap
# ============================================================

library(ggplot2)

# 1. Convert scaled matrix to long format
figure9_long <- as.data.frame(
  as.table(figure9_scaled_matrix),
  stringsAsFactors = FALSE
)

names(figure9_long) <- c(
  "gene_name",
  "comparison",
  "scaled_coefficient"
)

# 2. Preserve selection-rank order
figure9_long$gene_name <- factor(
  figure9_long$gene_name,
  levels = rev(rownames(figure9_scaled_matrix))
)

figure9_long$comparison <- factor(
  figure9_long$comparison,
  levels = colnames(figure9_scaled_matrix)
)

dim(figure9_long)

table(figure9_long$comparison)

anyNA(figure9_long)

all(
  figure9_long$scaled_coefficient >= -1 &
    figure9_long$scaled_coefficient <= 1
)

length(unique(figure9_long$gene_name))
# ============================================================
# S9.8e — Figure 9
# Step 4: Refined heatmap layout
# ============================================================

library(ggplot2)

# Short display labels
figure9_long$comparison_display <- factor(
  as.character(figure9_long$comparison),
  levels = c(
    "Age / 10 years",
    "Stage II vs I",
    "Stage III vs I",
    "Stage IV vs I"
  ),
  labels = c(
    "Age / 10y",
    "II vs I",
    "III vs I",
    "IV vs I"
  )
)

figure9_plot_refined <- ggplot(
  figure9_long,
  aes(
    x = comparison_display,
    y = gene_name,
    fill = scaled_coefficient
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  scale_fill_gradient2(
    low = "#426C9B",
    mid = "#FFFFFF",
    high = "#B65B68",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    name = "Within-contrast\nscaled coefficient"
  ) +
  scale_x_discrete(
    position = "top",
    expand = expansion(add = 0)
  ) +
  scale_y_discrete(
    expand = expansion(add = 0)
  ) +
  labs(
    title = "Integrated Age–Stage coefficient patterns",
    subtitle = "Top 30 stable shared genes | TCGA-BRCA",
    x = NULL,
    y = NULL,
    caption = paste0(
      "Adjusted model: Age per 10 years; Stage II, III and IV versus Stage I.\n",
      "Genes ordered by balanced adjusted-FDR ranking; no clustering.\n",
      "Coefficients scaled within each contrast; colors are not comparable across contrasts."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      color = "#252A33"
    ),
    plot.subtitle = element_text(
      size = 10,
      color = "#606873",
      margin = margin(b = 14)
    ),
    axis.text.x = element_text(
      size = 10,
      face = "bold",
      color = "#30343B",
      margin = margin(b = 8)
    ),
    axis.text.y = element_text(
      size = 9.5,
      color = "#30343B"
    ),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(
      size = 9,
      face = "bold"
    ),
    legend.text = element_text(size = 9),
    plot.caption = element_text(
      size = 8,
      color = "#69717D",
      hjust = 0,
      lineheight = 1.4,
      margin = margin(t = 15)
    ),
    plot.margin = margin(15, 25, 20, 15)
  )

print(figure9_plot_refined)

# ============================================================
# S9.8f — Step 1
# Inspect project working directory and output structure
# ============================================================

getwd()

list.files(
  path = ".",
  all.files = FALSE
)

# Confirm required objects are available
required_s9_objects <- c(
  "integrated_master_final",
  "figure7_plot",
  "figure8_plot",
  "figure9_plot_refined",
  "figure8_metrics",
  "figure9_top30",
  "figure9_coef_matrix",
  "figure9_scaled_matrix",
  "figure9_scale_factors"
)

stopifnot(all(vapply(required_s9_objects, function(x) exists(x, envir = analysis_env, inherits = FALSE), logical(1))))
data.frame(
  object = required_s9_objects,
  exists = vapply(
    required_s9_objects,
    exists,
    logical(1),
    inherits = FALSE
  )
)
# Check essential results in the Global Environment

objects_to_export <- c(
  "integrated_master_final",
  "figure9_ranked",
  "figure9_top30",
  "figure8_metrics",
  "figure9_scale_factors"
)

stopifnot(all(vapply(objects_to_export, function(x) exists(x, envir = analysis_env, inherits = FALSE), logical(1))))
data.frame(
  object = objects_to_export,
  available = vapply(
    objects_to_export,
    function(x) {
      exists(x, envir = analysis_env, inherits = FALSE)
    },
    logical(1)
  ),
  row.names = NULL
)
# ============================================================
# S9.8f — Step 6
# Export 418 stable shared genes with ranking metadata
# ============================================================

stopifnot(exists("figure9_ranked", inherits = FALSE))
dir.create(s9_export_dir, recursive = TRUE, showWarnings = FALSE)

ranked_csv <- file.path(
  s9_export_dir,
  "S9_stable_shared_genes_ranked.csv"
)

# Validate source object
stopifnot(
  nrow(figure9_ranked) == 418L,
  ncol(figure9_ranked) == 39L,
  !anyDuplicated(figure9_ranked$ensembl_gene_id),
  !file.exists(ranked_csv)
)

# Export
write.csv(
  figure9_ranked,
  file = ranked_csv,
  row.names = FALSE,
  na = ""
)

# Read back and validate
ranked_check <- read.csv(
  ranked_csv,
  check.names = FALSE
)

data.frame(
  check = c(
    "File exists",
    "Dimensions match",
    "Gene IDs and order match"
  ),
  passed = c(
    file.exists(ranked_csv),
    identical(dim(ranked_check), dim(figure9_ranked)),
    identical(
      ranked_check$ensembl_gene_id,
      figure9_ranked$ensembl_gene_id
    )
  )
)

dim(ranked_check)
# ============================================================
# S9.8g-A — Step 1
# Verify integrated evidence for Results narrative
# ============================================================

stopifnot(
  nrow(integrated_master_final) == 25551L,
  !anyDuplicated(integrated_master_final$ensembl_gene_id)
)

evidence_counts <- data.frame(
  analysis = c(
    "Primary Age",
    "Primary Stage",
    "Adjusted Age",
    "Adjusted Stage",
    "Primary Both",
    "Adjusted Both"
  ),
  significant_genes = c(
    sum(integrated_master_final$age_primary_significant),
    sum(integrated_master_final$stage_primary_significant),
    sum(integrated_master_final$age_adjusted_significant),
    sum(integrated_master_final$stage_adjusted_significant),
    sum(
      integrated_master_final$age_primary_significant &
        integrated_master_final$stage_primary_significant
    ),
    sum(
      integrated_master_final$age_adjusted_significant &
        integrated_master_final$stage_adjusted_significant
    )
  )
)

print(evidence_counts, row.names = FALSE)

cat(
  "\nPrimary association groups:\n"
)
print(table(integrated_master_final$association_group))

cat(
  "\nAdjusted association groups:\n"
)
print(table(integrated_master_final$adjusted_association_group))


# Explicit result exports; historical files are never reused as new outputs.
pipeline_event(ctx, "result_exports")
write.csv(age_results_annotated, file.path(output_dir, "data_processed/age_analysis/TCGA_BRCA_age_results_annotated.csv"), row.names = FALSE)
saveRDS(age_results_annotated, file.path(output_dir, "data_processed/age_analysis/TCGA_BRCA_age_results_annotated.rds"))
saveRDS(fit_age_eb, file.path(output_dir, "data_processed/age_analysis/TCGA_BRCA_age_limma_ebayes_fit.rds"))
saveRDS(gene_annotation, file.path(output_dir, "data_processed/age_analysis/TCGA_BRCA_gene_annotation_GENCODEv36.rds"))
write.csv(stage_results_annotated, file.path(output_dir, "data_processed/stage_analysis/TCGA_BRCA_stage_results_annotated.csv"), row.names = FALSE)
saveRDS(stage_results_annotated, file.path(output_dir, "data_processed/stage_analysis/TCGA_BRCA_stage_results_annotated.rds"))
saveRDS(fit_stage_eb, file.path(output_dir, "data_processed/stage_analysis/TCGA_BRCA_stage_limma_ebayes_fit.rds"))
saveRDS(stage_contrast_results, file.path(output_dir, "data_processed/stage_analysis/TCGA_BRCA_stage_pairwise_results.rds"))
write.csv(joint_age_results, file.path(output_dir, "data_processed/joint_analysis/TCGA_BRCA_joint_age_results.csv"), row.names = TRUE)
write.csv(joint_stage_omnibus, file.path(output_dir, "data_processed/joint_analysis/TCGA_BRCA_joint_stage_results.csv"), row.names = TRUE)
saveRDS(fit_joint_eb, file.path(output_dir, "data_processed/joint_analysis/TCGA_BRCA_joint_limma_ebayes_fit.rds"))
saveRDS(design_joint, file.path(output_dir, "data_processed/joint_analysis/TCGA_BRCA_joint_design.rds"))
saveRDS(joint_meta, file.path(output_dir, "data_processed/joint_analysis/TCGA_BRCA_joint_metadata.rds"))
saveRDS(voom_joint, file.path(output_dir, "data_processed/joint_analysis/TCGA_BRCA_joint_voom.rds"))
write.csv(integrated_master_final, file.path(output_dir, "data_processed/S9_integrated_results/S9_integrated_master.csv"), row.names = FALSE)
saveRDS(integrated_master_final, file.path(output_dir, "data_processed/S9_integrated_results/S9_integrated_master.rds"))
write.csv(figure8_metrics, file.path(output_dir, "data_processed/S9_integrated_results/S9_figure8_stability_metrics.csv"), row.names = FALSE)
write.csv(figure9_top30, file.path(output_dir, "data_processed/S9_integrated_results/S9_figure9_top30_genes.csv"), row.names = FALSE)
write.csv(figure9_coef_matrix, file.path(output_dir, "data_processed/S9_integrated_results/S9_figure9_adjusted_coefficients.csv"), row.names = TRUE)
write.csv(figure9_scaled_matrix, file.path(output_dir, "data_processed/S9_integrated_results/S9_figure9_scaled_coefficients.csv"), row.names = TRUE)
write.csv(data.frame(comparison = names(figure9_scale_factors), max_absolute_coefficient = unname(figure9_scale_factors)), file.path(output_dir, "data_processed/S9_integrated_results/S9_figure9_scaling_factors.csv"), row.names = FALSE)
ggsave(file.path(output_dir, "figures/S9_integrated_results/Figure_07_primary_association_categories.png"), figure7_plot, width = 9, height = 6.8, dpi = 300)
ggsave(file.path(output_dir, "figures/S9_integrated_results/Figure_08_coefficient_concordance.png"), figure8_plot, width = 10, height = 8, dpi = 300)
ggsave(file.path(output_dir, "figures/S9_integrated_results/Figure_09_stable_shared_top30_heatmap.png"), figure9_plot_refined, width = 10, height = 8, dpi = 300)
saveRDS(list(integrated_master_final = integrated_master_final, figure9_ranked = figure9_ranked, figure9_top30 = figure9_top30, figure8_metrics = figure8_metrics, figure9_coef_matrix = figure9_coef_matrix, figure9_scaled_matrix = figure9_scaled_matrix, figure9_scale_factors = figure9_scale_factors, figure7_plot = figure7_plot, figure8_plot = figure8_plot, figure9_plot_refined = figure9_plot_refined),
  file.path(output_dir, "data_processed/S9_integrated_results/recovery/S9_session_objects_v06.rds"))
invisible(list(output_dir = output_dir, evidence_counts = evidence_counts))
}
