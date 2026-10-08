suppressPackageStartupMessages({library(catboost);library(data.table);library(pROC);library(MASS);library(jsonlite)})
setDTthreads(1)
REF<-Sys.getenv('CIPDS_REFERENCE_DIR');SENS<-Sys.getenv('CIPDS_SENSITIVITY_DIR')
OUT<-file.path(SENS,'outputs');dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
THREADS<-23L;EPS<-1e-6;OUTER_K<-5L;BASE_SEED<-20260831L
for(expr in parse(file.path(REF,'scripts/14_train_overall_stacked_model.R'))){
 if(is.call(expr)&&identical(expr[[1]],as.name('<-'))&&is.call(expr[[3]])&&identical(expr[[3]][[1]],as.name('function')))eval(expr)
}
d<-fread(file.path(REF,'outputs/hospital_patient_model_ready_129.csv'))
fd<-fread(file.path(REF,'outputs/first_day_sensitivity/hospital_first_day_inputs.csv'))
features<-fread(file.path(REF,'outputs/nested_a_primary_final_features.csv'))
registry<-fread(file.path(REF,'outputs/nested_a_primary_calibration_and_model_registry.csv'))
fdreg<-fread(file.path(REF,'outputs/first_day_sensitivity/first_day_model_registry.csv'))
outcomes<-c('Outcome_NutriMetab','Outcome_TumorBurden','Outcome_TreatComp')
shorts<-c('NM','TB','TC')
stopifnot(nrow(d)==75248L,identical(d$patient_key,fd$patient_key))
missing_dir<-file.path(OUT,'missingness');first_dir<-file.path(OUT,'first_day_transfer')
dir.create(missing_dir,showWarnings=FALSE);dir.create(first_dir,showWarnings=FALSE)
summary_rows<-list();paired_rows<-list();missing_predictions<-list();calibration_rows<-list()
metrics<-function(y,p,component,analysis,n_features){
 ro<-roc(y,p,direction='<',quiet=TRUE);ci<-as.numeric(ci.auc(ro,method='delong'))
 data.table(component,analysis,n=length(y),positive=sum(y),prevalence=mean(y),n_features,auc=as.numeric(auc(ro)),auc_lower=ci[1],auc_upper=ci[3],average_precision=average_precision(y,p),brier=mean((y-p)^2))
}
for(oi in seq_along(outcomes)){
 outcome_name<-outcomes[oi];short<-shorts[oi];selected<-features[outcome==outcome_name,feature];reg<-registry[outcome==outcome_name]
 train<-d[split_calendar_entry=='train' & !is.na(get(outcome_name))];val<-d[split_calendar_entry=='validation' & !is.na(get(outcome_name))];test<-d[split_calendar_entry=='test' & !is.na(get(outcome_name))]
 missing_matrix<-function(z){m<-1.0*is.na(as.matrix(z[,..selected]));colnames(m)<-selected;m}
 fit<-fit_base(missing_matrix(train),train[[outcome_name]],reg$depth,reg$learning_rate,reg$l2_leaf_reg,2026100700L+oi)
 catboost.save_model(fit,file.path(missing_dir,paste0(short,'_missing_indicators.cbm')))
 raw_val<-predict_base(fit,missing_matrix(val));lp<-qlogis(pmin(pmax(raw_val,EPS),1-EPS));y<-val[[outcome_name]]
 cal<-glm(y~lp,family=binomial());stopifnot(all(is.finite(coef(cal))))
 raw_test<-predict_base(fit,missing_matrix(test));pm<-apply_platt(raw_test,coef(cal)[1],coef(cal)[2])
 primary<-catboost.load_model(file.path(REF,'outputs',reg$model_file))
 pp<-apply_platt(predict_base(primary,test[,..selected]),reg$intercept,reg$slope)
 yy<-test[[outcome_name]];complete<-rowSums(is.na(test[,..selected]))==0
 summary_rows[[length(summary_rows)+1L]]<-metrics(yy,pp,short,'Primary model; full test cohort',length(selected))
 summary_rows[[length(summary_rows)+1L]]<-metrics(yy,pm,short,'Missing indicators only; full test cohort',length(selected))
 summary_rows[[length(summary_rows)+1L]]<-metrics(yy[complete],pp[complete],short,'Primary model; all inputs observed',length(selected))
 rp<-roc(yy,pp,direction='<',quiet=TRUE);rm<-roc(yy,pm,direction='<',quiet=TRUE)
 delta<-as.numeric(auc(rp)-auc(rm));se<-sqrt(max(0,as.numeric(var(rp)+var(rm)-2*cov(rp,rm))))
 paired_rows[[oi]]<-data.table(component=short,n=length(yy),auc_primary=as.numeric(auc(rp)),auc_missing=as.numeric(auc(rm)),difference=delta,ci_lower=delta-qnorm(.975)*se,ci_upper=delta+qnorm(.975)*se,p_value=2*pnorm(-abs(delta/se)))
 calibration_rows[[oi]]<-data.table(component=short,validation_n=nrow(val),intercept=coef(cal)[1],slope=coef(cal)[2],features=paste(selected,collapse=';'))
 missing_predictions[[oi]]<-data.table(patient_key=test$patient_key,component=short,observed=yy,primary_probability=pp,missingness_probability=pm,observed_inputs=rowSums(!is.na(test[,..selected])))
 ft<-fd[match(test$patient_key,fd$patient_key)];fr<-fdreg[outcome==outcome_name]
 ff<-catboost.load_model(file.path(REF,'outputs/first_day_sensitivity',paste0(outcome_name,'_first_day.cbm')))
 fp<-apply_platt(predict_base(ff,ft[,..selected]),fr$intercept,fr$slope)
 gate<-rowSums(!is.na(ft[,..selected]))>0
 summary_rows[[length(summary_rows)+1L]]<-metrics(yy[gate],fp[gate],short,'First-day model; at least one input observed',length(selected))
 cat(short,'missingness diagnostic and complete-input evaluation done\n');flush.console()
}
paired<-rbindlist(paired_rows);paired[,p_holm:=p.adjust(p_value,'holm')]
fwrite(rbindlist(summary_rows),file.path(missing_dir,'hospital_diagnostic_performance.csv'))
fwrite(paired,file.path(missing_dir,'primary_vs_missingness_paired_auc.csv'))
fwrite(rbindlist(missing_predictions),file.path(missing_dir,'test_predictions.csv'))
fwrite(rbindlist(calibration_rows),file.path(missing_dir,'calibration_registry.csv'))

# Conditional cross-fitting: keep the already selected primary features and
# hyperparameters, as in the established first-day sensitivity analysis.
oof_list<-list();fold_checks<-list()
for(oi in seq_along(outcomes)){
 outcome_name<-outcomes[oi];short<-shorts[oi];selected<-features[outcome==outcome_name,feature]
 reg<-fdreg[outcome==outcome_name]
 train<-fd[split_calendar_entry=='train' & !is.na(get(outcome_name))];setorder(train,patient_key)
 fold_id<-make_stratified_folds(train[[outcome_name]],5L,BASE_SEED+oi)
 raw<-rep(NA_real_,nrow(train))
 for(k in 1:5){
  tr<-which(fold_id!=k);ho<-which(fold_id==k)
  stopifnot(length(intersect(train$patient_key[tr],train$patient_key[ho]))==0)
  fit<-fit_base(train[tr,..selected],train[[outcome_name]][tr],reg$depth,reg$learning_rate,reg$l2_leaf_reg,BASE_SEED+oi*10000L+k*1000L+900L)
  raw[ho]<-predict_base(fit,train[ho,..selected])
  fold_checks[[length(fold_checks)+1L]]<-data.table(component=short,fold=k,train_n=length(tr),holdout_n=length(ho),holdout_positives=sum(train[[outcome_name]][ho]),holdout_excluded_from_fit=TRUE)
  cat('First-day OOF',short,k,'/5 complete\n');flush.console()
 }
 stopifnot(!anyNA(raw))
 oo<-data.table(patient_key=train$patient_key,observed=train[[outcome_name]],probability=apply_platt(raw,reg$intercept,reg$slope))
 setnames(oo,c('observed','probability'),paste0(c('observed_','p_'),short));oof_list[[short]]<-oo
}
meta<-Reduce(function(a,b)merge(a,b,by='patient_key'),oof_list)
meta[,burden:=observed_NM+observed_TB+observed_TC]
scaling<-rbindlist(lapply(shorts,function(short){v<-qlogis(pmin(pmax(meta[[paste0('p_',short)]],EPS),1-EPS));data.table(component=short,mean=mean(v),sd=sd(v))}))
for(short in shorts){sr<-scaling[component==short];meta[[paste0('z_',short)]]<-(qlogis(pmin(pmax(meta[[paste0('p_',short)]],EPS),1-EPS))-sr$mean)/sr$sd}
meta[,burden_ordered:=ordered(burden,levels=0:3)]
model<-polr(burden_ordered~z_NM+z_TB+z_TC,data=meta,method='logistic',Hess=TRUE)
stopifnot(model$convergence==0L,nrow(meta)==49082L)
bundle<-list(model=model,component_order=shorts,logit_scaling=as.data.frame(scaling),phenoage_used=FALSE,nhanes_outcomes_used=FALSE,chosen_model_name='main_effects',training_population='Fivefold OOF predictions conditional on fixed primary feature sets and hyperparameters; first-day inputs')
saveRDS(bundle,file.path(first_dir,'overall_stacked_ordinal_model.rds'))
fwrite(meta[,.(patient_key,observed_NM,observed_TB,observed_TC,p_NM,p_TB,p_TC,burden)],file.path(first_dir,'overall_meta_training_oof.csv'))
fwrite(rbindlist(fold_checks),file.path(first_dir,'oof_fold_checks.csv'));fwrite(scaling,file.path(first_dir,'overall_logit_scaling.csv'))

mapped<-as.data.table(readRDS(file.path(REF,'outputs/nhanes_expanded_candidate_matrix.rds')))
score_variant<-function(variant,feature_table,reg_table,model_dir,ordinal_bundle){
 ans<-mapped[,.(SEQN,CYCLE)]
 for(oi in seq_along(outcomes)){
  outcome_name<-outcomes[oi];short<-shorts[oi];selected<-feature_table[outcome==outcome_name,feature];reg<-reg_table[outcome==outcome_name]
  fname<-if(variant=='first_day')paste0(outcome_name,'_first_day.cbm') else reg$model_file
  model<-catboost.load_model(file.path(model_dir,fname))
  pred<-apply_platt(predict_base(model,mapped[,..selected]),reg$intercept,reg$slope)
  nobs<-rowSums(!is.na(mapped[,..selected]));pred[nobs==0]<-NA_real_
  ans[[paste0(short,'_component')]]<-pred;ans[[paste0(short,'_input_observed_n')]]<-nobs
 }
 valid<-complete.cases(ans[,.(NM_component,TB_component,TC_component)])
 nd<-as.data.frame(ans[valid]);sc<-as.data.table(ordinal_bundle$logit_scaling)
 for(short in shorts){sr<-sc[component==short];nd[[paste0('z_',short)]]<-(qlogis(pmin(pmax(nd[[paste0(short,'_component')]],EPS),1-EPS))-sr$mean)/sr$sd}
 prob<-as.matrix(predict(ordinal_bundle$model,newdata=nd,type='probs'))[,as.character(0:3),drop=FALSE]
 ans[,Overall_expected_burden:=NA_real_];ans[valid,Overall_expected_burden:=as.numeric(prob%*%0:3)]
 stopifnot(nrow(ans)==44772L,uniqueN(ans$SEQN)==44772L,all(ans$Overall_expected_burden[valid]>=0 & ans$Overall_expected_burden[valid]<=3))
 fwrite(ans,file.path(OUT,paste0(variant,'_nhanes_scores.csv')))
}
score_variant('first_day',features,fdreg,file.path(REF,'outputs/first_day_sensitivity'),bundle)
V<-file.path(SENS,'variants/no_mpv/outputs')
nf<-fread(file.path(V,'nested_a_primary_final_features.csv'));nr<-fread(file.path(V,'nested_a_primary_calibration_and_model_registry.csv'))
stopifnot(!'lab_mpv'%in%nf$feature)
score_variant('no_mpv',nf,nr,V,readRDS(file.path(V,'overall_stacked_ordinal_model.rds')))
write_json(list(completed=as.character(Sys.time()),workers=23,first_day_oof_n=nrow(meta),first_day_feature_search_repeated=FALSE,first_day_hyperparameter_search_repeated=FALSE,first_day_models_reused=TRUE,no_mpv_selection_and_tuning_repeated=TRUE,component_observation_gate='At least one observed selected input',labels_changed=FALSE,primary_models_changed=FALSE),file.path(OUT,'hospital_sensitivity_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
cat('Hospital sensitivity diagnostics and both frozen NHANES score sets complete\n')
