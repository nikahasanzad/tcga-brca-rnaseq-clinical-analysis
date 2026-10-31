# Base-R checks of the released scientific expressions, output plan and inert entry point.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
release <- normalizePath(args[1], mustWork = TRUE)
source(file.path(release, "tests", "contract_tools.R"))
pipeline <- file.path(release, "R", "pipeline.R")
outer <- parse(pipeline, keep.source = FALSE)
stopifnot(length(outer) == 1L, identical(outer[[1]][[1]], as.name("<-")),
          identical(outer[[1]][[3]][[1]], as.name("function")))
statements <- pipeline_statements(pipeline)
get_assign <- function(name) Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                                    identical(x[[2]], as.name(name)), statements)
reference <- dget(file.path(release, "tests", "fixtures", "scientific_contract.R"))
for (name in names(reference$assignments)) {
  stopifnot(identical(get_assign(name), reference$assignments[[name]]))
}
for (name in names(reference$requests)) {
  request <- get_assign(name)[[1]][[3]][[4]]
  stopifnot(identical(as.list(request)[-1], reference$requests[[name]]))
}
lines <- readLines(pipeline, warn = FALSE)
get_block <- function(start, end) {
  a <- grep(start, lines, fixed = TRUE); b <- grep(end, lines, fixed = TRUE)
  stopifnot(length(a) == 1L, length(b) == 1L, a < b)
  parse(text = lines[a:(b - 1L)], keep.source = FALSE)
}
stopifnot(
  identical(get_block("# STEP 8 —", "# STEP 8.6 —"), reference$canonical_preprocessing),
  identical(get_block("# S9.5a — Joint cohort using canonical patient IDs", "# S9.5d — Joint Age + Stage inference"), reference$canonical_joint),
  identical(get_block("# STEP 2L — CANONICAL", "# STEP 2R — Final aliquot-level selection"), reference$technical_qc),
  identical(get_block("# S5.1 — CANONICAL", "# S5.2 —"), reference$contrasts)
)
cat("PASS: model formulas, fits, tests, preprocessing, contrasts, QC selection, ranking and final plots match frozen v05 expressions\n")
contract <- pipeline_contract(pipeline)
declared <- read.csv(file.path(release, "output_contract.csv"), stringsAsFactors = FALSE)
stopifnot(identical(contract, declared), nrow(contract) == 53L,
          !anyDuplicated(contract$relative_path))
for (name in c("read.delim", "read.csv", "readRDS", "write.csv", "saveRDS", "writeLines", "ggsave")) {
  definitions <- get_assign(name)
  stopifnot(length(definitions) == 1L,
            any(c("pipeline_read", "pipeline_publish") %in% all.names(definitions[[1]][[3]])))
}
symbols <- all.names(as.expression(statements))
stopifnot(!any(c("untar", "GDCquery", "GDCquery_clinic", "GDCdownload", "install.packages", ".GlobalEnv") %in% symbols),
          !any(grepl(paste0("/", "Users", "/"), lines, fixed = TRUE)))
calls <- function(name) Filter(function(x) is.call(x) && identical(x[[1]], as.name(name)), statements)
stopifnot(length(calls("pipeline_assert_fit")) == 3L,
          length(calls("pipeline_assert_results")) == 4L,
          length(calls("pipeline_assert_counts")) >= 1L)
cat("PASS: 53 unique output contracts, tracked read/write wrappers, guarded acquisition and critical validation handoffs\n")
cat("STRUCTURE_TESTS_COMPLETE; no scientific expressions were evaluated\n")
