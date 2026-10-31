args<-commandArgs(TRUE)
usage<-"Rscript --vanilla tools/restore_environment.R --library NEW_LIBRARY_PATH --execute"
if(!length(args)||"--help"%in%args||!"--execute"%in%args){cat(usage,"\n");quit(status=0)}
stopifnot(length(args)==3L,args[1]=="--library",args[3]=="--execute")
stopifnot(identical(as.character(getRversion()),"4.5.1"))
if(!requireNamespace("renv",quietly=TRUE))stop("Install the verified renv 1.3.1 packaging tool first; see environment/README.md.")
stopifnot(identical(as.character(packageVersion("renv")),"1.3.1"))
target<-normalizePath(args[2],mustWork=FALSE)
stopifnot(!file.exists(target),dir.exists(dirname(target)))
renv::restore(project=normalizePath("."),library=target,lockfile="renv.lock",
  clean=FALSE,prompt=FALSE)
