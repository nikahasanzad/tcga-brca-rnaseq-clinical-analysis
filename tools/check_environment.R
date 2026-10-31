args <- commandArgs(trailingOnly=TRUE)
if (length(args)) stop("No arguments or installation actions are supported.")
expected <- utils::read.csv("environment/direct_dependencies.csv", stringsAsFactors=FALSE)
cat("Recorded R: 4.5.1; installed R:", as.character(getRversion()), "\n")
ok <- identical(as.character(getRversion()), "4.5.1")
installed <- utils::installed.packages()
for (i in seq_len(nrow(expected))) {
  p <- expected$package[i]
  actual <- if (p %in% rownames(installed)) installed[p, "Version"] else NA_character_
  matched <- !is.na(actual) && identical(actual, expected$recorded_version[i])
  cat(p, "expected", expected$recorded_version[i], "installed", actual, if(matched) "MATCH" else "MISMATCH", "\n")
  ok <- ok && matched
}
cat("Metadata inspection only; no packages installed or scientific analysis run.\n")
if (!ok) quit(status=1)
