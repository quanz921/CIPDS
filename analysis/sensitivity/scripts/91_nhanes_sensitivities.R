REF<-Sys.getenv('CIPDS_REFERENCE_DIR');SENS<-Sys.getenv('CIPDS_SENSITIVITY_DIR')
source(file.path(REF,'scripts/27_supplement_common.R'))
suppressPackageStartupMessages({library(parallel);library(jsonlite);library(splines)})
setDTthreads(1);OUT<-file.path(SENS,'outputs');Q<-file.path(SENS,'qa')
x<-cipds_load_data(FALSE);original<-copy(x$dat)
make_design<-function(d)svydesign(ids=~survey_psu,strata=~survey_strata,weights=~pooled_mec_weight,data=d,nest=TRUE)
raw_names<-c(NM='NM_component',TB='TB_component',TC='TC_component',OVERALL='Overall_expected_burden')
variants<-list(primary=copy(original))
for(v in c('first_day','no_mpv')){
 d<-copy(original);s<-fread(file.path(OUT,paste0(v,'_nhanes_scores.csv')));idx<-match(d$SEQN,s$SEQN)
 stopifnot(!anyNA(idx),all(as.character(d$CYCLE)==as.character(s$CYCLE[idx])))
 for(raw in raw_names)d[[raw]]<-s[[raw]][idx]
 des<-make_design(d);agedes<-des[x$age60,]
 for(code in names(raw_names)){
  raw<-raw_names[[code]];mu<-as.numeric(coef(svymean(as.formula(paste0('~',raw)),agedes,na.rm=TRUE)));sd<-sqrt(as.numeric(svyvar(as.formula(paste0('~',raw)),agedes,na.rm=TRUE)))
  d[[CIPDS_SCORE_VARS[[code]]]]<-(d[[raw]]-mu)/sd
 }
 variants[[v]]<-d
}

# Reproduce the exact existing Q1-Q4 assignment, exposing the source sample,
# raw-unit and standardized cutpoints and observed counts for each cause.
common<-cipds_complete_mask(original,unname(CIPDS_SCORE_VARS));qdat<-original[common]
stopifnot(nrow(qdat)==9319L)
qcuts<-list();qcounts<-list();qr<-list()
qraw<-c(raw_names,PHENO='PhenoAge_acceleration')
for(code in names(CIPDS_SCORE_VARS)){
 z<-qdat[[CIPDS_SCORE_VARS[[code]]]];raw<-qdat[[qraw[[code]]]]
 cz<-cipds_weighted_quantile(z,qdat$pooled_mec_weight,c(.25,.5,.75));cr<-cipds_weighted_quantile(raw,qdat$pooled_mec_weight,c(.25,.5,.75))
 group<-cut(z,c(-Inf,cz,Inf),labels=FALSE,include.lowest=TRUE)
 raw_group<-cut(raw,c(-Inf,cr,Inf),labels=FALSE,include.lowest=TRUE)
 stopifnot(identical(group,raw_group))
 qcuts[[code]]<-data.table(score=code,reference_n=nrow(qdat),reference_deaths=sum(qdat$Death_AllCause),quantile=c(.25,.5,.75),cutpoint_raw=cr,cutpoint_standardized=cz,raw_variable=qraw[[code]],weight='Corrected pooled MEC weights',reference='Age >=60 common complete-case incremental cohort, all cancer-history groups combined')
 qr[[code]]<-data.table(SEQN=qdat$SEQN,score=code,quartile=group)
 for(cause in c('CANCER','CVD')){
  status<-qdat[[if(cause=='CANCER')'competing_cancer' else 'competing_cvd']];eligible<-!is.na(status)
  for(k in 1:4){keep<-group==k & eligible
   qcounts[[length(qcounts)+1L]]<-data.table(score=code,cause=cause,quartile=k,n=sum(keep),target_events=sum(status[keep]==1),competing_events=sum(status[keep]==2),censored=sum(status[keep]==0),weighted_percent=100*sum(qdat$pooled_mec_weight[keep])/sum(qdat$pooled_mec_weight[eligible]))
  }
 }
}
fwrite(rbindlist(qcuts),file.path(OUT,'quartile_cutpoints.csv'));fwrite(rbindlist(qcounts),file.path(OUT,'quartile_counts.csv'));fwrite(rbindlist(qr),file.path(OUT,'quartile_assignments.csv'))

# Compare mortality associations on the same participants for every variant.
assoc_mask<-Reduce('&',lapply(variants,function(d)cipds_complete_mask(d,unname(CIPDS_SCORE_VARS[names(raw_names)]))))
assoc_rows<-list()
for(v in names(variants)){
 d<-variants[[v]];des<-make_design(d)[assoc_mask,]
 for(code in names(raw_names)){
  score_var<-CIPDS_SCORE_VARS[[code]];fit<-svycoxph(cipds_make_formula(score_var),design=des)
  b<-coef(fit)[[score_var]];se<-sqrt(vcov(fit)[score_var,score_var])
  assoc_rows[[length(assoc_rows)+1L]]<-data.table(variant=v,score=code,n=sum(assoc_mask),events=sum(d$Death_AllCause[assoc_mask]),hr=exp(b),ci_lower=exp(b-qnorm(.975)*se),ci_upper=exp(b+qnorm(.975)*se),p_value=2*pnorm(-abs(b/se)))
 }
}
associations<-rbindlist(assoc_rows);associations[,p_holm:=p.adjust(p_value,'holm'),by=variant]
fwrite(associations,file.path(OUT,'variant_mortality_associations.csv'))

# The common reference cohort is held fixed before fitting any sensitivity model.
temporal_common<-Reduce('&',lapply(variants,function(d)cipds_complete_mask(d,unname(CIPDS_SCORE_VARS))))
train_mask<-temporal_common & original$CYCLE %in% c('2003-2004','2005-2006')
valid_mask<-temporal_common & original$CYCLE %in% c('2007-2008','2009-2010')
stopifnot(sum(train_mask)==2579L,sum(valid_mask)==3048L)
train_domain<-x$age60 & original$CYCLE %in% c('2003-2004','2005-2006')
pa_fit<-svyglm(PhenoAge~Age,design=x$design_full[train_domain,]);pa_coef<-coef(pa_fit)
resid<-original$PhenoAge-(pa_coef[1]+pa_coef[2]*original$Age)
transform_rows<-list();temporal_data<-list()
for(v in names(variants)){
 d<-copy(variants[[v]]);d[,pa_resid:=resid];des<-make_design(d)[train_domain,]
 for(code in names(CIPDS_SCORE_VARS)){
  raw<-if(code=='PHENO')'pa_resid' else raw_names[[code]]
  mu<-as.numeric(coef(svymean(as.formula(paste0('~',raw)),des,na.rm=TRUE)));sd<-sqrt(as.numeric(svyvar(as.formula(paste0('~',raw)),des,na.rm=TRUE)))
  d[[CIPDS_SCORE_VARS[[code]]]]<-(d[[raw]]-mu)/sd
  transform_rows[[length(transform_rows)+1L]]<-data.table(variant=v,score=code,training_mean=mu,training_sd=sd)
 }
 temporal_data[[v]]<-list(train=droplevels(copy(d[train_mask])),valid=droplevels(copy(d[valid_mask])),domain=d[train_domain])
}
fwrite(rbindlist(transform_rows),file.path(OUT,'temporal_training_transforms.csv'))
fwrite(data.table(term=names(pa_coef),coefficient=unname(pa_coef)),file.path(OUT,'temporal_phenoage_residualization.csv'))
grid<-1:10;predictions<-list();coefs<-list();baselines<-list();fit_objects<-list()
fit_one<-function(name,tr,va,formula){
 tr<-copy(tr);tr[,normalized_weight:=pooled_mec_weight/mean(pooled_mec_weight)]
 fit<-coxph(formula,data=tr,weights=normalized_weight,ties='efron',x=TRUE,y=TRUE,model=TRUE)
 stopifnot(all(is.finite(coef(fit))))
 predictions[[name]]<<-cipds_predict_risk(fit,va,grid);fit_objects[[name]]<<-fit
 coefs[[name]]<<-data.table(model=name,term=names(coef(fit)),coefficient=unname(coef(fit)))
 baselines[[name]]<<-cbind(data.table(model=name),as.data.table(basehaz(fit,centered=TRUE)))
}
linear_aug<-c('NM_z_age60','TB_z_age60','TC_z_age60','PHENO_z_age60')
fit_one('linear_base',temporal_data$primary$train,temporal_data$primary$valid,cipds_make_formula('PHENO_z_age60'))
for(v in names(temporal_data))fit_one(paste0(v,'_augmented'),temporal_data[[v]]$train,temporal_data[[v]]$valid,cipds_make_formula(linear_aug))

# Four knots define three natural-cubic basis terms per continuous covariate.
# All knot locations come exclusively from the early-period survey domain.
tr<-copy(temporal_data$primary$train);va<-copy(temporal_data$primary$valid);domain<-temporal_data$primary$domain
knot_rows<-list();spline_terms<-character()
for(variable in c('Age','PHENO_z_age60')){
 knots<-cipds_weighted_quantile(domain[[variable]],domain$pooled_mec_weight,c(.05,.35,.65,.95))
 stopifnot(length(unique(knots))==4L)
 bt<-ns(tr[[variable]],knots=knots[2:3],Boundary.knots=knots[c(1,4)],intercept=FALSE)
 bv<-predict(bt,newx=va[[variable]])
 cols<-paste0(if(variable=='Age')'age_ns_' else 'pheno_ns_',1:3)
 for(j in 1:3){tr[[cols[j]]]<-bt[,j];va[[cols[j]]]<-bv[,j]}
 spline_terms<-c(spline_terms,cols)
 knot_rows[[variable]]<-data.table(variable,percentile=c(.05,.35,.65,.95),knot=knots)
}
spline_clinical<-c(spline_terms,setdiff(CIPDS_CLINICAL_TERMS,'Age'))
form<-function(extra=character())as.formula(paste('Surv(Follow_Up_Years, Death_AllCause) ~',paste(c(extra,spline_clinical),collapse=' + ')))
fit_one('spline_base',tr,va,form());fit_one('spline_augmented',tr,va,form(c('NM_z_age60','TB_z_age60','TC_z_age60')))
fwrite(rbindlist(knot_rows),file.path(OUT,'nonlinear_training_knots.csv'))
fwrite(rbindlist(coefs),file.path(OUT,'temporal_frozen_coefficients.csv'));fwrite(rbindlist(baselines),file.path(OUT,'temporal_frozen_baselines.csv'))
saveRDS(fit_objects,file.path(OUT,'temporal_frozen_models.rds'))
valid<-temporal_data$primary$valid

# Reproduce the original pair on the same validation participants before using
# the shared evaluation routine for any new comparison.
reference<-readRDS(file.path(REF,'outputs/temporal_combination/validation_predictions.rds'))
stopifnot(identical(as.numeric(reference$SEQN),as.numeric(valid$SEQN)))
prediction_check<-rbindlist(lapply(c(linear_base='M2_CLINICAL_PHENO',primary_augmented='M5_CLINICAL_COMPONENTS_PHENO'),function(old)NULL))
for(pair in list(c('linear_base','M2_CLINICAL_PHENO'),c('primary_augmented','M5_CLINICAL_COMPONENTS_PHENO'))){
 err<-max(abs(predictions[[pair[1]]]$risk-reference$predictions[[pair[2]]]$risk));stopifnot(err<1e-9)
 cat('Primary frozen-risk reproduction:',pair[1],err,'\n')
}
registered<-parse(file.path(REF,'scripts/33_temporal_calibration.R'))
helper_names<-c('ipcw_binary_data','calibration_intercept_slope')
for(expr in registered)if(is.call(expr)&&identical(expr[[1]],as.name('<-'))&&is.symbol(expr[[2]])&&as.character(expr[[2]])%in%helper_names)eval(expr)
evaluate<-function(w){
 ans<-c();tt<-valid$Follow_Up_Years;ev<-valid$Death_AllCause
 for(m in names(predictions)){
  pred<-predictions[[m]];z<-c(uno_c10=cipds_weighted_uno(tt,ev,pred$lp,w,10))
  for(h in c(5,10)){
   z[paste0('auc',h)]<-cipds_weighted_td_auc(tt,ev,pred$lp,w,h)
   z[paste0('brier',h)]<-cipds_weighted_brier(tt,ev,pred$risk[,h],w,h)
   cal<-calibration_intercept_slope(tt,ev,pred$risk[,h],w,h)
   z[paste0('intercept',h)]<-cal['intercept'];z[paste0('slope',h)]<-cal['slope']
  }
  bg<-vapply(grid,function(h)cipds_weighted_brier(tt,ev,pred$risk[,h],w,h),numeric(1))
  z['ibs10']<-sum(diff(grid)*(head(bg,-1)+tail(bg,-1))/2)/9
  names(z)<-paste0(names(z),'__',m);ans<-c(ans,z)
 };ans
}
point<-evaluate(valid$pooled_mec_weight)
set.seed(2026091562L)
rep_full<-as.svrepdesign(x$design_full,type='bootstrap',replicates=1000,mse=TRUE);rep_sub<-rep_full[valid_mask,];rw<-weights(rep_sub,type='analysis')
cl<-makeCluster(23L)
clusterEvalQ(cl,{suppressPackageStartupMessages({library(data.table);library(survival)});setDTthreads(1);NULL})
clusterExport(cl,c('valid','predictions','grid','rw','evaluate','cipds_weighted_uno','cipds_weighted_km_censoring','cipds_weighted_td_auc','cipds_weighted_brier','ipcw_binary_data','calibration_intercept_slope'),envir=environment())
cat('Starting 1000 paired validation bootstrap replicates across six models with 23 workers\n');flush.console()
reps<-do.call(rbind,parLapplyLB(cl,1:1000,function(b)evaluate(rw[,b])))
stopCluster(cl);stopifnot(all(is.finite(reps)),identical(colnames(reps),names(point)))
df<-degf(x$design_full[valid_mask,]);crit<-qt(.975,df)
vv<-svrVar(reps,scale=rep_sub$scale,rscales=rep_sub$rscales,mse=rep_sub$mse,coef=point);se<-sqrt(diag(vv))
est<-rbindlist(lapply(seq_along(point),function(j){k<-strsplit(names(point)[j],'__',fixed=TRUE)[[1]];data.table(metric=k[1],model=k[2],n=nrow(valid),events=sum(valid$Death_AllCause),estimate=point[j],standard_error=se[j],ci_lower=point[j]-crit*se[j],ci_upper=point[j]+crit*se[j])}))
pairs<-list(primary=c('primary_augmented','linear_base'),first_day=c('first_day_augmented','linear_base'),no_mpv=c('no_mpv_augmented','linear_base'),nonlinear=c('spline_augmented','spline_base'))
comparisons<-rbindlist(lapply(names(pairs),function(v)rbindlist(lapply(c('uno_c10','auc5','auc10','brier5','brier10','ibs10'),function(metric){
 ka<-paste0(metric,'__',pairs[[v]][1]);kb<-paste0(metric,'__',pairs[[v]][2]);delta<-point[ka]-point[kb];rd<-reps[,ka]-reps[,kb]
 ds<-sqrt(as.numeric(svrVar(matrix(rd,ncol=1),scale=rep_sub$scale,rscales=rep_sub$rscales,mse=rep_sub$mse,coef=delta)));hb<-grepl('^(uno|auc)',metric)
 data.table(variant=v,metric,base_model=pairs[[v]][2],augmented_model=pairs[[v]][1],estimate_base=point[kb],estimate_augmented=point[ka],paired_difference=delta,ci_lower=delta-crit*ds,ci_upper=delta+crit*ds,benefit_difference=if(hb)delta else -delta,benefit_lower=if(hb)delta-crit*ds else -(delta+crit*ds),benefit_upper=if(hb)delta+crit*ds else -(delta-crit*ds),p_value=2*pt(-abs(delta/ds),df=df),higher_is_better=hb,n=nrow(valid),events=sum(valid$Death_AllCause))
}))))
comparisons[,p_holm:=p.adjust(p_value,'holm'),by=variant]
fwrite(est,file.path(OUT,'temporal_sensitivity_estimates.csv'));fwrite(comparisons,file.path(OUT,'temporal_sensitivity_comparisons.csv'))
saveRDS(list(point=point,replicates=reps,scale=rep_sub$scale,rscales=rep_sub$rscales,mse=rep_sub$mse,df=df),file.path(OUT,'temporal_paired_replicates.rds'))
saveRDS(list(SEQN=valid$SEQN,predictions=predictions,time=valid$Follow_Up_Years,event=valid$Death_AllCause,weights=valid$pooled_mec_weight),file.path(OUT,'temporal_validation_predictions.rds'))
cohorts<-data.table(analysis=c('Association','Temporal training','Temporal validation','Quartile reference'),n=c(sum(assoc_mask),sum(train_mask),sum(valid_mask),sum(common)),deaths=c(sum(original$Death_AllCause[assoc_mask]),sum(original$Death_AllCause[train_mask]),sum(original$Death_AllCause[valid_mask]),sum(original$Death_AllCause[common])))
fwrite(cohorts,file.path(OUT,'sensitivity_cohorts.csv'))
write_json(list(completed=as.character(Sys.time()),workers=23,bootstrap_replicates=1000,all_replicates_finite=TRUE,age60_source_n=sum(x$age60),full_survey_source_n=nrow(original),weights='Verified 4-year MEC weights for 1999-2002 and 2-year weights thereafter',temporal_training='2003-2006',temporal_validation='2007-2010',models_fixed_during_bootstrap=TRUE,holm='Six paired performance endpoints separately within each specification',spline='Natural cubic spline for age and PhenoAge Acceleration; four training-domain weighted knots at 5/35/65/95%; basis and residualization frozen',quartiles='Common 9319 complete cases; all cancer-history groups combined; cutpoints fixed across cancer and CVD analyses; cause availability filters applied after assignment'),file.path(OUT,'nhanes_sensitivity_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
print(comparisons[metric%in%c('uno_c10','ibs10')]);cat('NHANES sensitivity analyses complete\n')
