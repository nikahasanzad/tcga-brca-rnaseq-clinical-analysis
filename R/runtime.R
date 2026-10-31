# Runtime controls only; sourcing this file performs no I/O or analysis.
pipeline_require <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}
pipeline_inside <- function(path, root) {
  identical(path, root) || startsWith(path, paste0(root, .Platform$file.sep))
}
pipeline_md5 <- function(path) {
  pipeline_require(file.exists(path) && !dir.exists(path), paste("Missing file:", path))
  value <- unname(tools::md5sum(path))
  pipeline_require(length(value) == 1L && !is.na(value), paste("Cannot checksum:", path))
  value
}
pipeline_config <- function(input_dir, runs_dir, run_id, cache_dir = NULL,
                            allow_network = FALSE, seed = 42L) {
  pipeline_require(dir.exists(input_dir), "input_dir must already exist.")
  input_dir <- normalizePath(input_dir, mustWork = TRUE)
  pipeline_require(dir.exists(runs_dir), "runs_dir must already exist; choose a separate output parent.")
  runs_dir <- normalizePath(runs_dir, mustWork = TRUE)
  pipeline_require(!pipeline_inside(runs_dir, input_dir) && !pipeline_inside(input_dir, runs_dir),
                   "Input and run-output roots must be separate, non-nested directories.")
  pipeline_require(length(run_id) == 1L && grepl("^[A-Za-z0-9][A-Za-z0-9_-]{0,79}$", run_id),
                   "run_id must contain 1-80 letters/digits/underscores/hyphens; no path components.")
  pipeline_require(is.logical(allow_network) && length(allow_network) == 1L && !is.na(allow_network),
                   "allow_network must be TRUE or FALSE.")
  pipeline_require(length(seed) == 1L && is.finite(seed) && seed >= 0 &&
                     seed <= .Machine$integer.max && seed == as.integer(seed), "Invalid integer RNG seed.")
  if (is.null(cache_dir)) cache_dir <- file.path(input_dir, "gdc_cache")
  cache_dir <- normalizePath(cache_dir, mustWork = FALSE)
  pipeline_require(grepl("^/|^[A-Za-z]:", cache_dir), "cache_dir must be an absolute path.")
  list(input_dir = input_dir, runs_dir = runs_dir, run_id = run_id,
       cache_dir = cache_dir, allow_network = allow_network, seed = as.integer(seed),
       profile = "historical_TCGA_BRCA", checksum_algorithm = "MD5",
       formulas = c(age = "~ age_years", stage = "~ stage_group", joint = "~ age_10y + stage_group"),
       normalization = "TMM", filtering = "design-aware filterByExpr; package defaults recorded by version",
       primary_FDR = 0.05, stage_reference = "Stage I", age_joint_unit = "10 years",
       RNG_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
       qc_exclusions = c("TCGA-AC-A2QH", "TCGA-AC-A3OD", "TCGA-AC-A3QQ",
                         "TCGA-A7-A13G", "TCGA-A7-A26F", "TCGA-A7-A26I"))
}
pipeline_validate_config <- function(config) {
  keys <- c("input_dir", "runs_dir", "run_id", "cache_dir", "allow_network", "seed")
  pipeline_require(is.list(config) && all(keys %in% names(config)), "Invalid configuration schema.")
  expected <- do.call(pipeline_config, config[keys])
  pipeline_require(identical(config, expected), "Unsupported or altered settings; scientific settings must match the declared historical profile.")
}
pipeline_preflight <- function(config) {
  pipeline_validate_config(config)
  clinical <- file.path(config$input_dir, "data_raw", "clinical", "clinical.tsv")
  pipeline_require(file.exists(clinical), paste("Extracted clinical.tsv is required; no automatic untar:", clinical))
  for (name in c("rna_seq_technical_qc", "rna_seq_final")) {
    pipeline_require(dir.exists(file.path(config$input_dir, "data_raw", name)) || config$allow_network,
                     paste("Missing local count directory (network is disabled):", name))
  }
  if (!config$allow_network) for (name in c("query_exp", "biospecimen_raw")) {
    path <- file.path(config$cache_dir, paste0(name, ".rds"))
    pipeline_require(file.exists(path), paste("Offline GDC response cache required:", path,
      "Use a frozen request/value envelope or explicitly enable network for a new run."))
  }
  invisible(TRUE)
}
pipeline_new_context <- function(config, code_files, expected_outputs,
                                 required_packages = c("TCGAbiolinks", "edgeR", "limma", "ggplot2", "ggrepel", "scales")) {
  pipeline_validate_config(config)
  missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  pipeline_require(length(missing) == 0L, paste("Missing packages; install separately:", paste(missing, collapse = ", ")))
  pipeline_require(!anyDuplicated(basename(code_files)) && all(file.exists(code_files)), "Missing/ambiguous code files.")
  pipeline_require(length(expected_outputs) > 0L && !anyDuplicated(expected_outputs) &&
                     !any(grepl("^/|(^|/)\\.\\.(/|$)", expected_outputs)), "Invalid expected-output contract.")
  run_dir <- file.path(config$runs_dir, config$run_id)
  pipeline_require(!file.exists(run_dir) && !dir.exists(run_dir), paste("Refusing existing run:", run_dir))
  pipeline_require(dir.create(run_dir, recursive = FALSE), paste("Cannot reserve new run:", run_dir))
  ctx <- new.env(parent = emptyenv())
  ctx$config <- config; ctx$input_dir <- config$input_dir
  ctx$output_dir <- normalizePath(run_dir, mustWork = TRUE)
  ctx$stage <- "initialization"; ctx$inputs <- list(); ctx$outputs <- list(); ctx$assets <- list()
  ctx$events <- list(); ctx$requests <- list(); ctx$expected_outputs <- expected_outputs
  ctx$started_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  ctx$code <- lapply(code_files, function(p) list(path = normalizePath(p), md5 = pipeline_md5(p)))
  names(ctx$code) <- basename(code_files)
  ctx$packages <- setNames(lapply(required_packages, function(p) as.character(utils::packageVersion(p))), required_packages)
  ctx$network <- list(query = function(...) TCGAbiolinks::GDCquery(...),
                      clinic = function(...) TCGAbiolinks::GDCquery_clinic(...),
                      download = function(...) TCGAbiolinks::GDCdownload(...))
  pipeline_publish(ctx, "metadata/start.rds", function(p) base::saveRDS(
    list(config = config, code = ctx$code, packages = ctx$packages,
         R = R.version.string, started_utc = ctx$started_utc), p), record = FALSE)
  for (p in code_files) pipeline_publish(ctx, file.path("metadata", "code", basename(p)),
    function(dest) pipeline_require(file.copy(p, dest, overwrite = FALSE), "Cannot snapshot code."), record = FALSE)
  ctx
}
pipeline_event <- function(ctx, stage, detail = "") {
  ctx$stage <- stage
  message(format(Sys.time(), tz = "UTC", usetz = TRUE), " Stage: ", stage, if (nzchar(detail)) paste0(" — ", detail) else "")
  ctx$events[[length(ctx$events) + 1L]] <- list(time_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                                              stage = stage, detail = detail)
  invisible(NULL)
}
pipeline_target <- function(ctx, path) {
  pipeline_require(is.character(path) && length(path) == 1L && !is.na(path), "Invalid output path.")
  if (!grepl("^/|^[A-Za-z]:", path)) path <- file.path(ctx$output_dir, path)
  pipeline_require(!any(strsplit(gsub("\\\\", "/", path), "/", fixed = TRUE)[[1]] == ".."),
                   "Output traversal is forbidden.")
  # Validate the nearest existing ancestor before creating any directories.
  ancestor <- dirname(path)
  while (!dir.exists(ancestor)) {
    parent <- dirname(ancestor)
    pipeline_require(!identical(parent, ancestor), "Cannot resolve output parent.")
    ancestor <- parent
  }
  pipeline_require(pipeline_inside(normalizePath(ancestor), ctx$output_dir), "Output escapes run directory (possibly symlink).")
  link <- Sys.readlink(path)
  pipeline_require(!file.exists(path) && !dir.exists(path) && (is.na(link) || !nzchar(link)),
                   paste("Refusing existing output:", path))
  pipeline_require(dir.exists(dirname(path)) || dir.create(dirname(path), recursive = TRUE), "Cannot create output parent.")
  path <- file.path(normalizePath(dirname(path)), basename(path))
  pipeline_require(pipeline_inside(path, ctx$output_dir), "Output escapes run directory.")
  path
}
pipeline_publish <- function(ctx, path, writer, producer = "runtime metadata", settings = list(), record = TRUE) {
  target <- pipeline_target(ctx, path)
  # Preserve extension for ggsave's device detection. Hard-link commit cannot overwrite.
  temporary <- tempfile(".pending-", tmpdir = dirname(target), fileext = paste0(".", tools::file_ext(target)))
  on.exit(unlink(temporary), add = TRUE)
  writer(temporary)
  pipeline_require(file.exists(temporary) && file.info(temporary)$size > 0, "Writer produced no nonempty file.")
  pipeline_require(!file.exists(target) && isTRUE(file.link(temporary, target)),
                   paste("Atomic no-overwrite publication failed:", target))
  if (record) {
    row <- list(path = target, relative_path = substring(target, nchar(ctx$output_dir) + 2L),
                bytes = unname(file.info(target)$size), md5 = pipeline_md5(target),
                producer = producer, stage = ctx$stage, settings = settings,
                created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
    ctx$outputs[[row$relative_path]] <- row
    pipeline_publish(ctx, paste0(target, ".provenance.rds"), function(p) base::saveRDS(
      list(output = row, config = ctx$config, code = ctx$code, packages = ctx$packages,
           inputs_consumed_so_far = ctx$inputs, downloaded_inputs_consumed = ctx$assets,
           input_scope = "all registered inputs consumed before this output; may include upstream diagnostics"), p), record = FALSE)
  }
  invisible(target)
}
pipeline_track_input <- function(ctx, path) {
  pipeline_require(file.exists(path) && !dir.exists(path), paste("Missing input:", path))
  path <- normalizePath(path, mustWork = TRUE)
  roots <- c(ctx$input_dir, ctx$config$cache_dir, ctx$output_dir)
  pipeline_require(any(vapply(roots, function(root) pipeline_inside(path, root), logical(1))), "Read outside declared roots.")
  # Readbacks are output artifacts, not original inputs.
  if (pipeline_inside(path, ctx$output_dir)) return(path)
  info <- file.info(path)
  if (is.null(ctx$inputs[[path]])) {
    ctx$inputs[[path]] <- list(path = path, bytes = unname(info$size),
      mtime = as.numeric(info$mtime), md5 = pipeline_md5(path))
  } else pipeline_require(identical(ctx$inputs[[path]]$bytes, unname(info$size)) &&
      identical(ctx$inputs[[path]]$mtime, as.numeric(info$mtime)), paste("Input changed during run:", path))
  path
}
pipeline_read <- function(ctx, path, reader, ...) {
  path <- pipeline_track_input(ctx, path)
  value <- reader(path, ...)
  pipeline_track_input(ctx, path)
  value
}
pipeline_query <- function(ctx, key, request, kind) {
  path <- file.path(ctx$config$cache_dir, paste0(key, ".rds"))
  if (file.exists(path)) {
    envelope <- pipeline_read(ctx, path, base::readRDS)
    pipeline_require(is.list(envelope) && identical(envelope$schema_version, 1L) &&
      identical(envelope$request, request) && "value" %in% names(envelope), paste("Cache request/schema mismatch:", path))
    value <- envelope$value; origin <- "frozen cache"
  } else {
    pipeline_require(ctx$config$allow_network, paste("Network disabled; missing cache:", path))
    value <- do.call(ctx$network[[kind]], request); origin <- "live API"
    envelope <- list(schema_version = 1L, request = request, value = value)
    pipeline_publish(ctx, file.path("cache", paste0(key, ".rds")),
                     function(p) base::saveRDS(envelope, p), producer = paste("GDC", kind), settings = request)
  }
  ctx$requests[[key]] <- list(request = request, origin = origin, cache_input = if (origin == "frozen cache") path else NULL)
  value
}
pipeline_count_files <- function(ctx, query, results, folder, files_per_chunk) {
  pipeline_require(all(c("file_name", "file_id", "md5sum") %in% names(results)), "Count-file metadata schema changed.")
  pipeline_assert_ids(results$file_name, "count filenames")
  pipeline_require(all(basename(results$file_name) == results$file_name), "Invalid remote filenames.")
  pipeline_require(!anyNA(results$md5sum) && all(grepl("^[0-9a-fA-F]{32}$", results$md5sum)), "Missing/invalid expected count MD5.")
  local_dir <- file.path(ctx$input_dir, "data_raw", folder)
  local <- if (dir.exists(local_dir)) list.files(local_dir, recursive = TRUE, full.names = TRUE) else character()
  relevant <- local[basename(local) %in% results$file_name]
  pipeline_require(!anyDuplicated(basename(relevant)), "Ambiguous duplicate local count filenames.")
  missing <- setdiff(results$file_name, basename(relevant))
  if (length(missing)) {
    pipeline_require(ctx$config$allow_network, paste("Network disabled; missing count files:", paste(missing, collapse = ", ")))
    destination <- file.path(ctx$output_dir, "downloads", folder)
    pipeline_require(!dir.exists(destination), "Download staging directory already exists.")
    dir.create(destination, recursive = TRUE)
    reduced <- query
    reduced$results[[1]] <- query$results[[1]][query$results[[1]]$file_name %in% missing, , drop = FALSE]
    ctx$network$download(query = reduced, method = "api", files.per.chunk = files_per_chunk, directory = destination)
    added <- list.files(destination, recursive = TRUE, full.names = TRUE)
    added <- added[basename(added) %in% missing]
    pipeline_require(!anyDuplicated(basename(added)) && setequal(basename(added), missing), "Downloaded count-file set incomplete/ambiguous.")
    relevant <- c(relevant, added)
    pipeline_event(ctx, ctx$stage, paste("Downloaded", length(added), "missing files to isolated staging"))
  }
  pipeline_require(setequal(basename(relevant), results$file_name), "Count-file set mismatch.")
  relevant <- sort(relevant)
  for (path in relevant) {
    pipeline_track_input(ctx, path)
    expected <- results$md5sum[match(basename(path), results$file_name)]
    pipeline_require(identical(tolower(pipeline_md5(path)), tolower(expected)), paste("Count checksum mismatch:", path))
    if (pipeline_inside(normalizePath(path), ctx$output_dir)) {
      ctx$assets[[normalizePath(path)]] <- list(path = normalizePath(path),
        bytes = unname(file.info(path)$size), md5 = pipeline_md5(path),
        producer = "GDCdownload", expected_md5 = expected,
        file_id = results$file_id[match(basename(path), results$file_name)], stage = ctx$stage)
    }
  }
  relevant
}
pipeline_assert_ids <- function(ids, label = "IDs") {
  pipeline_require(is.character(ids) && length(ids) > 0L && !anyNA(ids) &&
    all(nzchar(ids)) && !anyDuplicated(ids), paste("Missing/empty/duplicate", label))
}
pipeline_assert_counts <- function(x) {
  pipeline_require(is.matrix(x) && is.numeric(x) && nrow(x) > 0 && ncol(x) > 0 &&
    all(is.finite(x)) && all(x >= 0) && all(x == round(x)) && all(colSums(x) > 0), "Invalid raw count matrix.")
  pipeline_assert_ids(rownames(x), "gene IDs"); pipeline_assert_ids(colnames(x), "patient IDs")
}
pipeline_assert_fit <- function(fit, voom, design) {
  pipeline_require(identical(rownames(fit$coefficients), rownames(voom$E)) &&
    identical(colnames(fit$coefficients), colnames(design)) &&
    identical(rownames(design), colnames(voom$E)) &&
    all(is.finite(fit$coefficients)) && all(is.finite(fit$sigma)) && all(fit$sigma >= 0) &&
    all(fit$df.residual > 0) && all(fit$df.residual == nrow(design) - ncol(design)), "Fit/design/gene alignment or residual degrees of freedom failed.")
}
pipeline_assert_results <- function(x) {
  pipeline_assert_ids(rownames(x), "result gene IDs")
  pipeline_require(all(c("P.Value", "adj.P.Val") %in% names(x)), "Result schema missing P/FDR.")
  p <- x$P.Value; q <- x$adj.P.Val
  pipeline_require(all(is.finite(p)) && all(is.finite(q)) && all(p >= 0 & p <= 1) &&
    all(q >= 0 & q <= 1) && isTRUE(all.equal(p.adjust(p, "BH"), q, tolerance = 1e-12)), "Invalid P/FDR or BH correction mismatch.")
}
pipeline_finalize <- function(ctx, status, error = NULL) {
  pipeline_require(status %in% c("COMPLETE", "FAILED"), "Invalid run status.")
  if (status == "COMPLETE") {
    pipeline_require(all(ctx$expected_outputs %in% names(ctx$outputs)), "Required outputs missing; run cannot be COMPLETE.")
    for (row in ctx$inputs) pipeline_require(identical(pipeline_md5(row$path), row$md5), paste("Input bytes changed:", row$path))
    for (row in ctx$outputs) pipeline_require(identical(pipeline_md5(row$path), row$md5), paste("Output bytes changed:", row$path))
    for (row in ctx$assets) pipeline_require(identical(pipeline_md5(row$path), row$md5), paste("Downloaded input bytes changed:", row$path))
    for (row in ctx$code) pipeline_require(identical(pipeline_md5(row$path), row$md5), "Code changed during run.")
  }
  final_dir <- file.path("metadata", tolower(status))
  pipeline_publish(ctx, file.path(final_dir, "sessionInfo.txt"), function(p) writeLines(capture.output(sessionInfo()), p), record = FALSE)
  manifest <- list(schema_version = 1L, status = status, error = error,
    started_utc = ctx$started_utc, finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    config = ctx$config, code = ctx$code, packages = ctx$packages, inputs = ctx$inputs,
    outputs = ctx$outputs, downloaded_inputs = ctx$assets, requests = ctx$requests, events = ctx$events,
    expected_outputs = ctx$expected_outputs,
    provenance_scope = "Per-output sidecars list all inputs consumed so far; final manifest verifies all original inputs and output checksums.")
  pipeline_publish(ctx, file.path(final_dir, "manifest.rds"), function(p) base::saveRDS(manifest, p), record = FALSE)
  for (field in c("inputs", "outputs")) {
    items <- manifest[[field]]
    table <- if (length(items)) do.call(rbind, lapply(items, function(x) {
      data.frame(path = x$path, bytes = x$bytes, md5 = x$md5, stringsAsFactors = FALSE)
    })) else data.frame(path = character(), bytes = numeric(), md5 = character())
    pipeline_publish(ctx, file.path(final_dir, paste0(field, ".csv")),
                     function(p) utils::write.csv(table, p, row.names = FALSE), record = FALSE)
  }
  pipeline_publish(ctx, "STATUS.txt", function(p) writeLines(status, p), record = FALSE)
  invisible(manifest)
}

# ggplot objects created inside a function retain that function's full environment.
# Replace only the analysis-frame reference in a copied object; retain plot data,
# expressions and package namespaces. The on-disk RDS then needs no custom reader.
pipeline_detach_analysis <- function(object, analysis_env, bindings = character()) {
  stopifnot(is.environment(analysis_env))
  portable_env <- new.env(parent = baseenv())
  pipeline_require(is.character(bindings) && !anyDuplicated(bindings), "Invalid portable plot bindings.")
  for (name in bindings) {
    pipeline_require(exists(name, envir = analysis_env, inherits = FALSE), paste("Missing portable plot binding:", name))
    value <- get(name, envir = analysis_env, inherits = FALSE)
    pipeline_require(is.atomic(value) && as.numeric(object.size(value)) < 1024^2,
                     paste("Portable plot bindings must be small atomic display metadata:", name))
    assign(name, value, envir = portable_env)
  }
  marker <- "TCGA_BRCA_DETACHED_ANALYSIS_FRAME_v1"
  bytes <- serialize(object, NULL, refhook = function(x) {
    if (is.environment(x) && identical(x, analysis_env)) marker else NULL
  })
  unserialize(bytes, refhook = function(key) {
    pipeline_require(identical(key, marker), "Unknown serialization reference.")
    portable_env
  })
}

# Native macOS PDF supports Unicode without requiring an XQuartz installation.
# Other systems use Cairo; successful portability testing is required separately.
pipeline_pdf_device <- function(filename, width, height, bg = "white", ...) {
  if (identical(Sys.info()[["sysname"]], "Darwin")) {
    grDevices::quartz(type = "pdf", file = filename, width = width, height = height, bg = bg, ...)
  } else {
    grDevices::cairo_pdf(filename = filename, width = width, height = height, bg = bg, ...)
  }
}
