# Isolated runtime tests only. All fixture files live under this test's temp directory.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
release <- normalizePath(args[1], mustWork = TRUE)
source(file.path(release, "R", "runtime.R"))
source(file.path(release, "run.R"))
temporary <- tempfile("phase32-runtime-")
dir.create(temporary)
input <- file.path(temporary, "input"); runs <- file.path(temporary, "runs")
dir.create(input); dir.create(runs)
dir.create(file.path(input, "data_raw", "clinical"), recursive = TRUE)
writeLines("patient_id\nfixture", file.path(input, "data_raw", "clinical", "clinical.tsv"))
for (d in c("rna_seq_technical_qc", "rna_seq_final")) dir.create(file.path(input, "data_raw", d))
cache <- file.path(input, "gdc_cache"); dir.create(cache)
for (name in c("query_exp", "biospecimen_raw")) saveRDS(list(schema_version = 1L, request = list(), value = list()), file.path(cache, paste0(name, ".rds")))
expect_error <- function(expr) stopifnot(inherits(try(force(expr), silent = TRUE), "try-error"))
context <- function(id, expected = "results.rds", online = FALSE) {
  pipeline_new_context(pipeline_config(input, runs, id, allow_network = online),
    file.path(release, "R", "runtime.R"), expected, required_packages = character())
}

ctx <- context("publication")
altered_config <- pipeline_config(input, runs, "altered-settings")
altered_config$primary_FDR <- .5
expect_error(pipeline_validate_config(altered_config))
expect_error(pipeline_config(input, input, "nested"))
expect_error(pipeline_config(input, runs, "../escape"))
data_file <- file.path(input, "original.txt"); writeLines("before", data_file)
invisible(pipeline_read(ctx, data_file, readLines))
pipeline_event(ctx, "fixture")
pipeline_publish(ctx, "results.rds", function(p) saveRDS(list(value = 3), p),
                 producer = "fixture_object", settings = list(scale = "test"))
digest <- pipeline_md5(file.path(ctx$output_dir, "results.rds"))
expect_error(pipeline_publish(ctx, "results.rds", function(p) stop("writer must not run")))
expect_error(pipeline_publish(ctx, file.path(input, "overwrite.txt"), function(p) writeLines("bad", p)))
expect_error(pipeline_publish(ctx, "../escape.txt", function(p) writeLines("bad", p)))
stopifnot(!file.exists(file.path(input, "overwrite.txt")), identical(digest, pipeline_md5(file.path(ctx$output_dir, "results.rds"))))
if (file.symlink(input, file.path(ctx$output_dir, "symlink"))) {
  expect_error(pipeline_publish(ctx, "symlink/escape.txt", function(p) writeLines("bad", p)))
}
manifest <- pipeline_finalize(ctx, "COMPLETE")
sidecar <- readRDS(file.path(ctx$output_dir, "results.rds.provenance.rds"))
stopifnot(identical(sidecar$output$producer, "fixture_object"),
          identical(sidecar$output$settings$scale, "test"),
          length(sidecar$inputs_consumed_so_far) == 1L,
          identical(manifest$status, "COMPLETE"),
          identical(readLines(file.path(ctx$output_dir, "STATUS.txt")), "COMPLETE"))
expect_error(context("publication"))
cat("PASS: isolated output roots, no overwrite/traversal/symlink escape, sidecars and complete manifest\n")

race <- context("race")
expect_error(pipeline_publish(race, "results.rds", function(p) {
  saveRDS("temporary", p)
  writeLines("competing file", file.path(race$output_dir, "results.rds"))
}))
stopifnot(identical(readLines(file.path(race$output_dir, "results.rds")), "competing file"))
cat("PASS: a destination created during writing is not overwritten\n")

changed <- context("changed-input")
invisible(pipeline_read(changed, data_file, readLines))
pipeline_publish(changed, "results.rds", function(p) saveRDS(1, p))
old_time <- file.info(data_file)$mtime
writeLines("after!", data_file); Sys.setFileTime(data_file, old_time)
expect_error(pipeline_finalize(changed, "COMPLETE"))
pipeline_finalize(changed, "FAILED", "fixture input mutation")
stopifnot(identical(readLines(file.path(changed$output_dir, "STATUS.txt")), "FAILED"))
tampered <- context("changed-output")
pipeline_publish(tampered, "results.rds", function(p) saveRDS(1, p))
saveRDS(2, file.path(tampered$output_dir, "results.rds"))
expect_error(pipeline_finalize(tampered, "COMPLETE"))
cat("PASS: changed input bytes and tampered outputs cannot be marked complete\n")

offline <- context("offline")
request <- list(project = "fixture")
saveRDS(list(schema_version = 1L, request = request, value = list(ok = TRUE)), file.path(cache, "fixture.rds"))
offline$network$query <- function(...) stop("Network must never be called")
stopifnot(isTRUE(pipeline_query(offline, "fixture", request, "query")$ok))
expect_error(pipeline_query(offline, "fixture", list(project = "different"), "query"))
expect_error(pipeline_query(offline, "missing", request, "query"))
online <- context("online-fixture", online = TRUE)
online$network$query <- function(...) list(mock = TRUE)
stopifnot(isTRUE(pipeline_query(online, "new", request, "query")$mock),
          file.exists(file.path(online$output_dir, "cache", "new.rds")),
          !file.exists(file.path(cache, "new.rds")))
cat("PASS: frozen cache request checks; offline blocks network; opt-in cache writes stay inside the new run\n")

local_file <- file.path(input, "data_raw", "rna_seq_technical_qc", "a.tsv")
writeLines("counts-A", local_file)
staged_fixture <- file.path(temporary, "b.tsv"); writeLines("counts-B", staged_fixture)
results <- data.frame(file_name = c("a.tsv", "b.tsv"), file_id = c("A", "B"),
                      md5sum = c(pipeline_md5(local_file), pipeline_md5(staged_fixture)))
query <- list(results = list(results))
expect_error(pipeline_count_files(offline, query, results, "rna_seq_technical_qc", 10L))
online$network$download <- function(query, method, files.per.chunk, directory) {
  stopifnot(identical(query$results[[1]]$file_name, "b.tsv"))
  file.copy(staged_fixture, file.path(directory, "b.tsv"))
}
files <- pipeline_count_files(online, query, results, "rna_seq_technical_qc", 10L)
stopifnot(setequal(basename(files), results$file_name), length(online$assets) == 1L,
          identical(pipeline_md5(local_file), results$md5sum[1]))
wrong <- results[1, , drop = FALSE]; wrong$md5sum <- paste(rep("0", 32), collapse = "")
expect_error(pipeline_count_files(offline, query, wrong, "rna_seq_technical_qc", 10L))
cat("PASS: exact local file selection, count MD5 checks, missing-only downloads and downloaded-input provenance\n")

expect_error(pipeline_assert_ids(c("a", "a")))
expect_error(pipeline_assert_counts(matrix(c(1, -1), nrow = 1)))
design <- cbind("(Intercept)" = 1, age = 1:4)
rownames(design) <- paste0("p", 1:4)
voom <- list(E = matrix(1:12, nrow = 3, dimnames = list(paste0("g", 1:3), rownames(design))))
fit <- list(coefficients = matrix(1:6, nrow = 3,
  dimnames = list(rownames(voom$E), colnames(design))), sigma = rep(1, 3), df.residual = rep(2, 3))
pipeline_assert_fit(fit, voom, design)
fit$df.residual[1] <- 0
expect_error(pipeline_assert_fit(fit, voom, design))
p <- c(.01, .2, .5)
r <- data.frame(P.Value = p, adj.P.Val = p.adjust(p, "BH"), row.names = c("g1", "g2", "g3"))
pipeline_assert_results(r); r$adj.P.Val[1] <- .9
expect_error(pipeline_assert_results(r))
cat("PASS: invalid IDs, counts and inconsistent BH-FDR stop\n")

set.seed(902)
seed <- .Random.seed
success_config <- pipeline_config(input, runs, "runner-fixture")
success <- run_pipeline(success_config, release, required_packages = character(),
  analysis_function = function(ctx) {
    pipeline_track_input(ctx, local_file)
    for (path in ctx$expected_outputs) pipeline_publish(ctx, path,
      function(p) writeLines("synthetic runtime fixture; not a scientific result", p), producer = "synthetic runtime fixture")
    invisible(NULL)
  })
stopifnot(identical(.Random.seed, seed), length(success$manifest$outputs) == 53L,
          identical(success$manifest$status, "COMPLETE"))
failure_config <- pipeline_config(input, runs, "runner-failure")
expect_error(run_pipeline(failure_config, release, required_packages = character(),
  analysis_function = function(ctx) stop("intentional fixture failure")))
failed <- readRDS(file.path(runs, "runner-failure", "metadata", "failed", "manifest.rds"))
stopifnot(identical(failed$status, "FAILED"), grepl("intentional fixture failure", failed$error))
missing_config <- pipeline_config(input, runs, "runner-missing-output")
expect_error(run_pipeline(missing_config, release, required_packages = character(),
  analysis_function = function(ctx) invisible(NULL)))
stopifnot(identical(readLines(file.path(runs, "runner-missing-output", "STATUS.txt")), "FAILED"))
cat("PASS: injected synthetic runner completion/failure, missing-output refusal, console capture and RNG restoration\n")

before <- list.files(temporary, recursive = TRUE, all.files = TRUE)
inert <- new.env(parent = globalenv())
sys.source(file.path(release, "R", "runtime.R"), inert)
sys.source(file.path(release, "R", "pipeline.R"), inert)
sys.source(file.path(release, "run.R"), inert)
invisible(capture.output(inert$pipeline_cli(c("--input-dir", "/not-used", "--run-id", "no-execute"))))
stopifnot(identical(before, list.files(temporary, recursive = TRUE, all.files = TRUE)))
cat("PASS: sourcing and CLI without --execute do not start analysis or create run files\n")
unlink(temporary, recursive = TRUE)
cat("RUNTIME_TESTS_COMPLETE; all data and completion markers above are disposable synthetic fixtures\n")
