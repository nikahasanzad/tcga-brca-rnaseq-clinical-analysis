args<-commandArgs(TRUE);root<-normalizePath(if(length(args))args[1]else".",mustWork=TRUE)
rscript<-file.path(R.home("bin"),"Rscript")
commands<-list(c("tools/check_integrity.R","."),c("tests/test_structure.R","."),
  c("tests/test_runtime.R","."),c("tests/test_optional_qc.R","."),
  c("tests/test_public_tables.R","."),c("run.R","--help"),
  c("tools/prepare_archived_cache.R","--help"))
old<-setwd(root);on.exit(setwd(old))
for(command in commands){
  cat("CHECK",command[1],"\n")
  status<-system2(rscript,c("--vanilla",command))
  if(status!=0L)stop("Offline release check failed: ",command[1])
}
cat("OFFLINE_RELEASE_TESTS_COMPLETE; no TCGA pipeline executed\n")
