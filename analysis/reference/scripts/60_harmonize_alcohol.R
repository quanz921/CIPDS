suppressPackageStartupMessages({library(data.table); library(jsonlite)})
root <- Sys.getenv('CIPDS_PACKAGE_DIR', 'cipds_revision_20260915')
e <- new.env(); load(file.path(root,'inputs/nhanes_phase7_batch1.RData'),envir=e)
d <- as.data.table(e$nhanes_valid); old <- copy(d)
stopifnot(nrow(d)==44772L, uniqueN(d$SEQN)==44772L)
spec <- data.table(CYCLE=c('1999-2000','2001-2002'),year=c(1999,2001),
  field=c('ALQ100','ALD100'),file=c('ALQ_1999.xpt','ALQ_B_2001.xpt'))
d[, Alcohol_12Drinks_AnyYear := as.numeric(ALQ101)]
d[, Alcohol_Source_Field := 'ALQ101']
for(i in seq_len(nrow(spec))) {
  s <- spec[i]
  raw <- fread(file.path(root,'inputs',paste0('alcohol_',s$year,'_original.csv')))
  stopifnot(uniqueN(raw$SEQN)==nrow(raw),s$field %in% names(raw))
  idx <- match(d$SEQN,raw$SEQN)
  v <- as.numeric(raw[[s$field]][idx])
  in_cycle <- d$CYCLE==s$CYCLE
  stopifnot(all(is.na(v[!in_cycle])))
  d[[s$field]] <- v
  d[in_cycle, Alcohol_12Drinks_AnyYear := v[in_cycle]]
  d[in_cycle, Alcohol_Source_Field := s$field]
}
d[, Alcohol_Status := fifelse(Alcohol_12Drinks_AnyYear==1,'Drinker',
  fifelse(Alcohol_12Drinks_AnyYear==2,'Non-drinker',NA_character_))]
d[, Alcohol_F := factor(Alcohol_Status,levels=levels(old$Alcohol_F))]
unchanged <- setdiff(names(old),c('Alcohol_F','Alcohol_Status'))
stopifnot(all(vapply(unchanged,function(v)identical(old[[v]],d[[v]]),logical(1))))
later <- !d$CYCLE %in% spec$CYCLE
stopifnot(identical(old$Alcohol_F[later],d$Alcohol_F[later]),identical(old$WTMEC2YR,d$WTMEC2YR))
d[, old_missing := is.na(old$Alcohol_F)]
mapping <- d[Age>=60,.(n=.N,raw_field=unique(Alcohol_Source_Field),
  valid_yes=sum(Alcohol_12Drinks_AnyYear==1,na.rm=TRUE),
  valid_no=sum(Alcohol_12Drinks_AnyYear==2,na.rm=TRUE),
  missing_before=sum(old_missing),missing_after=sum(is.na(Alcohol_F)),
  missing_pct_before=100*mean(old_missing),missing_pct_after=100*mean(is.na(Alcohol_F)),
  weighted_missing_pct_before=100*weighted.mean(old_missing,WTMEC2YR),
  weighted_missing_pct_after=100*weighted.mean(is.na(Alcohol_F),WTMEC2YR)),by=CYCLE]
mapping[, harmonized_field := 'Alcohol_12Drinks_AnyYear']
mapping[, recode := '1=yes; 2=no; 7/9/missing=missing']
fwrite(mapping,file.path(root,'outputs/alcohol_mapping_by_cycle.csv'),bom=TRUE)
changes <- d[,.(SEQN,CYCLE,Age,raw_field=Alcohol_Source_Field,
  raw_response=Alcohol_12Drinks_AnyYear,old_missing,new_status=as.character(Alcohol_F))]
fwrite(changes,file.path(root,'qa/alcohol_harmonization_record_audit.csv'))
d[,old_missing:=NULL]
nhanes_valid <- as.data.frame(d)
save(nhanes_valid,file=file.path(root,'inputs/nhanes_phase7_batch1.RData'),compress='gzip')
write_json(list(source_n=nrow(d),unchanged_existing_columns=length(unchanged),
  all_other_existing_columns_identical=TRUE,later_cycle_alcohol_identical=TRUE,
  weights_unchanged=TRUE,weights_scope='WTMEC2YR/9 retained at user request',
  alcohol_definition='Self-reported consumption of at least 12 drinks in any one year; cycle-specific questionnaire serving definitions',
  raw_fields=c('1999-2000: ALQ100','2001-2002: ALD100','2003-2016: ALQ101')),
  file.path(root,'qa/alcohol_harmonization_qa.json'),auto_unbox=TRUE,pretty=TRUE)
print(mapping)
