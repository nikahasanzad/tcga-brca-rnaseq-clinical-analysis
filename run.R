# Source defines functions only. CLI needs --execute to start the analysis.
run_pipeline <- function(config, release_dir, analysis_function = NULL,
                          required_packages = c("TCGAbiolinks", "edgeR", "limma", "ggplot2", "ggrepel", "scales")) {
  release_dir <- normalizePath(release_dir, mustWork = TRUE)
  env <- new.env(parent = globalenv())
  sys.source(file.path(release_dir, "R", "runtime.R"), envir = env)
  sys.source(file.path(release_dir, "R", "pipeline.R"), envir = env)
  env$pipeline_preflight(config)
  contract_file <- file.path(release_dir, "output_contract.csv")
  expected <- utils::read.csv(contract_file, stringsAsFactors = FALSE)$relative_path
  code_files <- c(file.path(release_dir, "run.R"), file.path(release_dir, "R", c("runtime.R", "pipeline.R")), contract_file)
  ctx <- env$pipeline_new_context(config, code_files, expected, required_packages)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(config$seed, kind = config$RNG_kind[1], normal.kind = config$RNG_kind[2], sample.kind = config$RNG_kind[3])
  grDevices::pdf(file = NULL)
  device <- grDevices::dev.cur()
  on.exit(if (device %in% grDevices::dev.list()) grDevices::dev.off(device), add = TRUE)
  # Keep default plotting away from the user's working directory.
  if (is.null(analysis_function)) analysis_function <- env$tcga_brca_analysis
  console <- character()
  connection <- textConnection("console", "w", local = TRUE)
  sink(connection)
  sink_open <- TRUE
  close_console <- function() {
    if (sink_open) { sink(); close(connection); sink_open <<- FALSE }
  }
  on.exit(close_console(), add = TRUE)
  tryCatch({
    value <- withCallingHandlers(analysis_function(ctx), warning = function(w) {
      env$pipeline_event(ctx, ctx$stage, paste("WARNING:", conditionMessage(w)))
    })
    close_console()
    env$pipeline_publish(ctx, "metadata/console.log", function(p) writeLines(c("Run console", console), p), record = FALSE)
    manifest <- env$pipeline_finalize(ctx, "COMPLETE")
    invisible(list(value = value, manifest = manifest, run_dir = ctx$output_dir))
  }, error = function(e) {
    close_console()
    if (!file.exists(file.path(ctx$output_dir, "metadata", "console.log"))) {
      try(env$pipeline_publish(ctx, "metadata/console.log", function(p) writeLines(c("Run console", console), p), record = FALSE), silent = TRUE)
    }
    env$pipeline_event(ctx, ctx$stage, paste("ERROR:", conditionMessage(e)))
    try(env$pipeline_finalize(ctx, "FAILED", conditionMessage(e)), silent = TRUE)
    stop(e)
  })
}

pipeline_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- paste(
    "No analysis runs without --execute.",
    "Rscript --vanilla run.R --input-dir /read-only/project --runs-dir /separate/existing/runs --run-id brca-001 --execute",
    "Optional: --cache-dir /frozen/gdc_cache --seed 42 --allow-network",
    "Offline is the default. No automatic installation or archive extraction.", sep = "\n")
  if (!length(args) || "--help" %in% args) { cat(usage, "\n"); return(invisible(NULL)) }
  options <- list(); flags <- c("--execute", "--allow-network")
  values <- c("--input-dir", "--runs-dir", "--run-id", "--cache-dir", "--seed")
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (!key %in% c(flags, values)) stop("Unknown argument: ", key, call. = FALSE)
    if (!is.null(options[[key]])) stop("Duplicate argument: ", key, call. = FALSE)
    if (key %in% flags) { options[[key]] <- TRUE; i <- i + 1L }
    else {
      if (i == length(args) || startsWith(args[i + 1L], "--")) stop("Missing value: ", key, call. = FALSE)
      options[[key]] <- args[i + 1L]; i <- i + 2L
    }
  }
  if (!isTRUE(options[["--execute"]])) { cat(usage, "\n"); return(invisible(NULL)) }
  for (key in c("--input-dir", "--runs-dir", "--run-id")) if (is.null(options[[key]])) stop("Required: ", key, call. = FALSE)
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) != 1L) stop("Run this entry point through Rscript.", call. = FALSE)
  release_dir <- dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
  controls <- new.env(parent = globalenv())
  sys.source(file.path(release_dir, "R", "runtime.R"), controls)
  seed <- if (is.null(options[["--seed"]])) 42 else suppressWarnings(as.numeric(options[["--seed"]]))
  config <- controls$pipeline_config(options[["--input-dir"]], options[["--runs-dir"]],
    options[["--run-id"]], cache_dir = options[["--cache-dir"]],
    allow_network = isTRUE(options[["--allow-network"]]), seed = seed)
  result <- run_pipeline(config, release_dir)
  cat("COMPLETE:", result$run_dir, "\n")
  invisible(result)
}

if (sys.nframe() == 0L) pipeline_cli()
