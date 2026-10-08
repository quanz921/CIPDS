source(file.path(Sys.getenv('CIPDS_PACKAGE_DIR','cipds_revision_20260915'),'scripts/27_supplement_common.R'))
suppressPackageStartupMessages({library(parallel);library(jsonlite)})
setDTthreads(1)
out<-file.path(CIPDS_ROOT,'outputs/temporal_combination');dir.create(out,recursive=TRUE,showWarnings=FALSE)
x<-cipds_load_data(FALSE);dat<-x$dat
common<-cipds_complete_mask(dat,unname(CIPDS_SCORE_VARS))
train_mask<-common & dat$CYCLE %in% c('2003-2004','2005-2006')
valid_mask<-common & dat$CYCLE %in% c('2007-2008','2009-2010')
train<-droplevels(copy(dat[train_mask]));valid<-droplevels(copy(dat[valid_mask]))
train_original<-copy(train);valid_original<-copy(valid)
train_domain<-x$age60 & dat$CYCLE %in% c('2003-2004','2005-2006')
design_train_domain<-x$design_full[train_domain,]
pheno_age_fit<-svyglm(PhenoAge~Age,design=design_train_domain)
pa_coef<-coef(pheno_age_fit)
pa_resid<-dat$PhenoAge-(pa_coef[1]+pa_coef[2]*dat$Age)
design_train_domain<-update(design_train_domain,pa_resid=pa_resid[train_domain])
scale_rows<-list()
for(code in c('NM','TB','TC','OVERALL','PHENO')) {
 raw_name<-switch(code,NM='NM_component',TB='TB_component',TC='TC_component',OVERALL='Overall_expected_burden',PHENO='pa_resid')
 mu<-as.numeric(coef(svymean(as.formula(paste0('~',raw_name)),design_train_domain,na.rm=TRUE)))
 sigma<-sqrt(as.numeric(svyvar(as.formula(paste0('~',raw_name)),design_train_domain,na.rm=TRUE)))
 stopifnot(is.finite(sigma),sigma>0)
 if(code=='PHENO') {
  train[[CIPDS_SCORE_VARS[[code]]]]<-(pa_resid[train_mask]-mu)/sigma
  valid[[CIPDS_SCORE_VARS[[code]]]]<-(pa_resid[valid_mask]-mu)/sigma
 } else {
  train[[CIPDS_SCORE_VARS[[code]]]]<-(train[[raw_name]]-mu)/sigma
  valid[[CIPDS_SCORE_VARS[[code]]]]<-(valid[[raw_name]]-mu)/sigma
 }
 scale_rows[[code]]<-data.table(code,raw_name,training_mean=mu,training_sd=sigma)
}
fwrite(rbindlist(scale_rows),file.path(out,'training_transform_registry.csv'))
fwrite(data.table(term=names(pa_coef),coefficient=unname(pa_coef)),file.path(out,'training_phenoage_age_regression.csv'))
models<-list(M2_CLINICAL_PHENO='PHENO_z_age60',
 M5_CLINICAL_COMPONENTS_PHENO=c('NM_z_age60','TB_z_age60','TC_z_age60','PHENO_z_age60'))
labels<-c(M2_CLINICAL_PHENO='Clinical base + PhenoAge Acceleration',
 M5_CLINICAL_COMPONENTS_PHENO='Clinical base + NM + TB + TC + PhenoAge Acceleration')
grid<-1:10;predictions<-list();identity<-list();coefficients<-list();baselines<-list()
for(m in names(models)) {
 train[,normalized_weight:=pooled_mec_weight/mean(pooled_mec_weight)]
 fit<-coxph(cipds_make_formula(models[[m]]),data=train,weights=normalized_weight,ties='efron',x=TRUE,y=TRUE,model=TRUE)
 predictions[[m]]<-cipds_predict_risk(fit,valid,grid)
 train_original[,normalized_weight:=pooled_mec_weight/mean(pooled_mec_weight)]
 oldfit<-coxph(cipds_make_formula(models[[m]]),data=train_original,weights=normalized_weight,ties='efron',x=TRUE,y=TRUE,model=TRUE)
 oldpred<-cipds_predict_risk(oldfit,valid_original,grid)
 err<-max(abs(predictions[[m]]$risk-oldpred$risk))
 stopifnot(err<1e-8)
 identity[[m]]<-data.table(model=m,max_absolute_risk_difference_from_prior_preprocessing=err)
 coefficients[[m]]<-data.table(model=m,term=names(coef(fit)),coefficient=unname(coef(fit)))
 baselines[[m]]<-cbind(data.table(model=m),as.data.table(basehaz(fit,centered=TRUE)))
}
fwrite(rbindlist(identity),file.path(out,'preprocessing_invariance_check.csv'))
fwrite(rbindlist(coefficients),file.path(out,'frozen_coefficients.csv'))
fwrite(rbindlist(baselines),file.path(out,'frozen_baseline_cumulative_hazard.csv'))
# Reuse the registered calibration estimators, without executing its analysis.
registered<-parse(file.path(CIPDS_ROOT,'scripts/33_temporal_calibration.R'))
helper_names<-c('ipcw_binary_data','calibration_intercept_slope')
for(expr in registered)if(is.call(expr)&&identical(expr[[1]],as.name('<-'))&&is.symbol(expr[[2]])&&as.character(expr[[2]])%in%helper_names)eval(expr)
evaluate<-function(w) {
 ans<-c();tt<-valid$Follow_Up_Years;ev<-valid$Death_AllCause
 for(m in names(models)) {
  pred<-predictions[[m]]
  z<-c(uno_c10=cipds_weighted_uno(tt,ev,pred$lp,w,10))
  for(h in c(5,10)) {
   z[paste0('auc',h)]<-cipds_weighted_td_auc(tt,ev,pred$lp,w,h)
   z[paste0('brier',h)]<-cipds_weighted_brier(tt,ev,pred$risk[,h],w,h)
   cal<-calibration_intercept_slope(tt,ev,pred$risk[,h],w,h)
   z[paste0('intercept',h)]<-cal['intercept'];z[paste0('slope',h)]<-cal['slope']
  }
  bg<-vapply(grid,function(h)cipds_weighted_brier(tt,ev,pred$risk[,h],w,h),numeric(1))
  z['ibs10']<-sum(diff(grid)*(head(bg,-1)+tail(bg,-1))/2)/9
  names(z)<-paste0(names(z),'__',m);ans<-c(ans,z)
 }
 ans
}
point<-evaluate(valid$pooled_mec_weight)
set.seed(2026091562L)
rep_full<-as.svrepdesign(x$design_full,type='bootstrap',replicates=1000,mse=TRUE)
rep_sub<-rep_full[valid_mask,];rw<-weights(rep_sub,type='analysis')
cl<-makeCluster(CIPDS_THREADS);on.exit(try(stopCluster(cl),silent=TRUE),add=TRUE)
clusterEvalQ(cl,{suppressPackageStartupMessages({library(data.table);library(survival)});setDTthreads(1);NULL})
clusterExport(cl,c('valid','models','predictions','grid','rw','evaluate','cipds_weighted_uno',
 'cipds_weighted_km_censoring','cipds_weighted_td_auc','cipds_weighted_brier',
 'ipcw_binary_data','calibration_intercept_slope'),envir=environment())
reps<-do.call(rbind,parLapplyLB(cl,1:1000,function(b)evaluate(rw[,b])))
stopCluster(cl);cl<-NULL
stopifnot(all(is.finite(reps)),identical(colnames(reps),names(point)))
df<-degf(x$design_full[valid_mask,]);crit<-qt(.975,df)
vv<-svrVar(reps,scale=rep_sub$scale,rscales=rep_sub$rscales,mse=rep_sub$mse,coef=point)
se<-sqrt(diag(vv))
est<-rbindlist(lapply(seq_along(point),function(j){k<-strsplit(names(point)[j],'__',fixed=TRUE)[[1]];
 data.table(metric=k[1],model=k[2],model_label=labels[[k[2]]],n=nrow(valid),events=sum(valid$Death_AllCause),
 estimate=point[j],standard_error=se[j],ci_lower=point[j]-crit*se[j],ci_upper=point[j]+crit*se[j])}))
comparisons<-rbindlist(lapply(c('uno_c10','auc5','auc10','brier5','brier10','ibs10'),function(metric){
 ka<-paste0(metric,'__M5_CLINICAL_COMPONENTS_PHENO');kb<-paste0(metric,'__M2_CLINICAL_PHENO')
 delta<-point[ka]-point[kb];rd<-reps[,ka]-reps[,kb]
 ds<-sqrt(as.numeric(svrVar(matrix(rd,ncol=1),scale=rep_sub$scale,rscales=rep_sub$rscales,mse=rep_sub$mse,coef=delta)))
 hb<-grepl('^(uno|auc)',metric)
 data.table(metric,estimate_base=point[kb],estimate_augmented=point[ka],paired_difference=delta,
  ci_lower=delta-crit*ds,ci_upper=delta+crit*ds,standard_error=ds,higher_is_better=hb,
  benefit_difference=if(hb)delta else -delta,p_value=2*pt(-abs(delta/ds),df=df),n=nrow(valid),events=sum(valid$Death_AllCause))
}))
comparisons[,p_holm:=p.adjust(p_value,'holm')]
fwrite(est,file.path(out,'temporal_combination_estimates.csv'),bom=TRUE)
fwrite(comparisons,file.path(out,'temporal_combination_paired_comparisons.csv'),bom=TRUE)
saveRDS(list(point=point,replicates=reps,scale=rep_sub$scale,rscales=rep_sub$rscales,mse=rep_sub$mse,df=df),file.path(out,'paired_replicates.rds'))
saveRDS(list(SEQN=valid$SEQN,predictions=predictions,time=valid$Follow_Up_Years,event=valid$Death_AllCause,weights=valid$pooled_mec_weight),file.path(out,'validation_predictions.rds'))
fwrite(data.table(role=c('Training','Temporal validation'),cycles=c('2003-2006','2007-2010'),
 n=c(nrow(train),nrow(valid)),deaths=c(sum(train$Death_AllCause),sum(valid$Death_AllCause)),
 maximum_follow_up=c(max(train$Follow_Up_Years),max(valid$Follow_Up_Years))),file.path(out,'cohort_registry.csv'))
write_json(list(training='2003-2006',validation='2007-2010',train_n=nrow(train),validation_n=nrow(valid),
 preprocessing='Training-period survey-domain score means/SDs and PhenoAge-on-age regression; frozen for validation',
 model='Unpenalized sampling-weighted Cox; coefficients and baseline hazard trained in early cycles',
 comparison='Clinical + PhenoAge Acceleration + NM/TB/TC versus clinical + PhenoAge Acceleration',
 uncertainty='1000 paired validation survey-bootstrap replicate weights; training models held fixed',
 multiplicity='Holm across six paired performance endpoints',primary_endpoint='10-year Uno C-index',
 no_15_year_extrapolation=TRUE,all_replicates_finite=TRUE,post_2002_weights_unchanged=TRUE),
 file.path(out,'analysis_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
print(comparisons)
