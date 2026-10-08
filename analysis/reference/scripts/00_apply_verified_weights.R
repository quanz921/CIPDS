suppressPackageStartupMessages({library(data.table);library(haven);library(jsonlite)})
setDTthreads(1)
root <- Sys.getenv('CIPDS_PACKAGE_DIR')
stopifnot(nzchar(root))
source(file.path(root,'scripts/weight_contract.R'))
input <- file.path(root,'inputs/nhanes_phase7_batch1.RData')
env <- new.env()
objects <- load(input,envir=env)
original <- as.data.frame(get('nhanes_valid',envir=env))
stopifnot(nrow(original)==44772L,!anyDuplicated(original$SEQN))
official <- rbindlist(lapply(c('1999-2000','2001-2002'),function(cycle){
  filename <- paste0('DEMO_',gsub('-','_',cycle),'.XPT')
  d <- as.data.table(read_xpt(file.path(root,'inputs/cdc_weight_audit_20260908',filename)))
  stopifnot(all(c('SEQN','WTMEC2YR','WTMEC4YR')%in%names(d)),!anyDuplicated(d$SEQN))
  data.table(SEQN=as.numeric(d$SEQN),CYCLE=cycle,official_WTMEC2YR=as.numeric(d$WTMEC2YR),official_WTMEC4YR=as.numeric(d$WTMEC4YR))
}))
stopifnot(!anyDuplicated(official$SEQN))
d <- copy(as.data.table(original))
early <- as.character(d$CYCLE)%in%c('1999-2000','2001-2002')
stopifnot(sum(early)==8983L)
idx <- match(d$SEQN[early],official$SEQN)
stopifnot(!anyNA(idx),all(as.character(d$CYCLE[early])==official$CYCLE[idx]))
stopifnot(all(abs(d$WTMEC2YR[early]-official$official_WTMEC2YR[idx])<1e-7))
d[,WTMEC4YR:=NA_real_]
d[which(early),WTMEC4YR:=official$official_WTMEC4YR[idx]]
d[,pooled_mec_weight:=cipds_pooled_mec_weights(d)]
old <- d$WTMEC2YR/9
stopifnot(all(d$pooled_mec_weight[!early]==old[!early]),all(d$WTMEC4YR[early]>0))
stopifnot(identical(as.data.frame(d[,names(original),with=FALSE]),original))
audit <- d[,.(n=.N,changed=sum(abs(pooled_mec_weight-WTMEC2YR/9)>1e-9),
              original_sum=sum(WTMEC2YR/9),corrected_sum=sum(pooled_mec_weight)),by=CYCLE]
fwrite(audit,file.path(root,'qa/weight_change_by_cycle.csv'))
fwrite(d[,.(SEQN,CYCLE,WTMEC2YR,WTMEC4YR,pooled_mec_weight)],file.path(root,'qa/verified_weight_registry.csv'))
replacement <- as.data.frame(d)
if (is.data.table(get('nhanes_valid',envir=env))) replacement <- d
assign('nhanes_valid',replacement,envir=env)
save(list=objects,envir=env,file=input,compress=TRUE)
result <- list(status='PASS',source_n=nrow(d),early_n=sum(early),official_matches=length(idx),
  added_columns=c('WTMEC4YR','pooled_mec_weight'),original_columns_unchanged=TRUE,
  original_WTMEC2YR_preserved=TRUE,post_2002_weights_unchanged=TRUE,
  formula=CIPDS_WEIGHT_DESCRIPTION,
  official_urls=c('https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/DEMO.xpt',
                  'https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/DEMO_B.xpt'))
write_json(result,file.path(root,'qa/weight_join_validation.json'),pretty=TRUE,auto_unbox=TRUE)
print(audit)
cat('Weight join and original-field preservation checks PASSED\n')
