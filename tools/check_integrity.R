args<-commandArgs(TRUE);root<-if(length(args))args[1]else"."
manifest<-readLines(file.path(root,"SHA256SUMS"),warn=FALSE)
stopifnot(length(manifest)>0L,all(grepl("^[0-9a-f]{64}  .+$",manifest)))
expected<-substr(manifest,1,64);paths<-substring(manifest,67)
stopifnot(!anyDuplicated(paths),!any(grepl("^/|(^|/)\\.\\.(/|$)",paths)))
exe<-Sys.which("sha256sum");extra<-character()
if(!nzchar(exe)){exe<-Sys.which("shasum");extra<-c("-a","256")}
stopifnot(nzchar(exe))
for(i in seq_along(paths)){
  f<-file.path(root,paths[i]);stopifnot(file.exists(f),!dir.exists(f))
  actual<-system2(exe,c(extra,shQuote(f)),stdout=TRUE,stderr=TRUE)
  stopifnot(is.null(attr(actual,"status")),length(actual)>=1L,identical(substr(actual[1],1,64),expected[i]))
}
cat("PASS: SHA256 for",length(paths),"manifest entries\n")
