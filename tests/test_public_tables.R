args<-commandArgs(TRUE);root<-if(length(args))args[1]else"."
spec<-data.frame(file=c("TCGA_BRCA_age_results_annotated.csv","TCGA_BRCA_stage_results_annotated.csv",
  "TCGA_BRCA_joint_age_results.csv","TCGA_BRCA_joint_stage_results.csv"),
  genes=c(25551,31157,31285,31285),significant=c(9063,1285,9413,1201))
for(i in seq_len(nrow(spec))){
  x<-read.csv(file.path(root,"tables","baseline",spec$file[i]),check.names=FALSE)
  stopifnot(nrow(x)==spec$genes[i],sum(x$adj.P.Val<0.05)==spec$significant[i])
}
x<-read.csv(file.path(root,"validation","V4_candidate_status_public.csv"),check.names=FALSE)
stopifnot(nrow(x)==2120L,!any(c("V2_worst_removed_patient_id","baseline_evidence_path","source_V2")%in%names(x)))
n<-table(x$classification)
stopifnot(n[["supported across the specified checks"]]==745L,
  n[["sensitivity-qualified"]]==103L,n[["insufficient information / not estimable"]]==1272L)
a<-read.csv(file.path(root,"validation","V4_claim_evidence_assessment.csv"))
stopifnot(sum(a$evidence_status=="ESTIMABLE_BASELINE_NOT_SIGNIFICANT")==1272L)
cat("PASS: archived table counts and V4 classifications; no P/FDR recalculation\n")
