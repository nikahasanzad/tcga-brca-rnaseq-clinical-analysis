# AST inspection only; never evaluate scientific statements or call writers.
pipeline_statements <- function(path) {
  ast <- parse(path, keep.source = FALSE)
  entry <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                    identical(x[[2]], as.name("tcga_brca_analysis")), as.list(ast))
  stopifnot(length(entry) == 1L)
  as.list(entry[[1]][[3]][[3]])[-1L]
}
pipeline_contract <- function(path) {
  statements <- pipeline_statements(path)
  root <- "/run"
  env <- list2env(list(output_dir = root, processed_dir = file.path(root, "data_processed"),
    age_figures_dir = file.path(root, "figures", "age_analysis"),
    s9_export_dir = file.path(root, "data_processed", "S9_integrated_results"),
    ranked_csv = file.path(root, "data_processed", "S9_integrated_results", "S9_stable_shared_genes_ranked.csv")),
    parent = baseenv())
  rows <- list()
  for (x in statements) {
    if (!is.call(x) || !is.symbol(x[[1]])) next
    fn <- as.character(x[[1]])
    if (!fn %in% c("write.csv", "saveRDS", "writeLines", "ggsave")) next
    args <- as.list(x)[-1]; n <- names(args)
    key <- if (fn == "ggsave") "filename" else if (fn == "writeLines") "con" else "file"
    pos <- if (!is.null(n) && key %in% n) match(key, n) else if (fn == "ggsave") 1L else 2L
    destination <- eval(args[[pos]], env)
    objkey <- if (fn == "ggsave") "plot" else if (fn == "writeLines") "text" else "x"
    objpos <- if (!is.null(n) && objkey %in% n) match(objkey, n) else if (fn == "ggsave") 2L else 1L
    stopifnot(startsWith(destination, "/run/"))
    rows[[length(rows) + 1L]] <- data.frame(function_name = fn,
      relative_path = substring(destination, nchar("/run/") + 1L),
      producer = paste(deparse(args[[objpos]]), collapse = " "), stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  stopifnot(!anyDuplicated(out$relative_path))
  out
}
