source(file.path(Sys.getenv("CIPDS_PACKAGE_DIR"),"scripts/27_supplement_common.R"))
x<-cipds_load_data(FALSE); d<-x$dat
masks<-list(source=x$age60,overall=x$age60 & !is.na(d$OVERALL_z_age60),
 components=cipds_complete_mask(d,c("NM_z_age60","TB_z_age60","TC_z_age60")),
 paired=x$age60 & complete.cases(d[,.(NM_z_age60,TB_z_age60,TC_z_age60,OVERALL_z_age60,PHENO_z_age60)]),
 common=cipds_complete_mask(d,unname(CIPDS_SCORE_VARS)))
r<-rbindlist(lapply(names(masks),function(k)data.table(cohort=k,n=sum(masks[[k]]),events=sum(d$Death_AllCause[masks[[k]]]))))
fwrite(r,file.path(CIPDS_ROOT,"qa/cohort_registry.csv"));print(r)
fwrite(d[,.(SEQN,CYCLE,Age,PHENO_z_age60,PhenoAge_acceleration,Diabetes,Hypertension,CVD,Follow_Up_Years)],file.path(CIPDS_ROOT,"qa/analysis_input_registry.csv"))
