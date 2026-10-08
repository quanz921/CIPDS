REF<-Sys.getenv('CIPDS_REFERENCE_DIR');SENS<-Sys.getenv('CIPDS_SENSITIVITY_DIR')
source(file.path(REF,'scripts/27_supplement_common.R'))
suppressPackageStartupMessages(library(jsonlite));setDTthreads(1)
OUT<-file.path(SENS,'outputs');V<-file.path(SENS,'variants/no_mpv/outputs');Q<-file.path(SENS,'qa')
checks<-list()
check<-function(name,pass,detail=''){
 checks[[length(checks)+1L]]<<-data.table(check=name,passed=isTRUE(pass),detail=as.character(detail))
 if(!isTRUE(pass))cat('FAILED',name,detail,'\n')
}
near<-function(a,b,tol=1e-9)length(a)==length(b)&&all(is.finite(a))&&all(is.finite(b))&&max(abs(a-b))<tol
x<-cipds_load_data(FALSE);d<-x$dat
check('Source survey sample',nrow(d)==44772L&&sum(x$age60)==15048L)
for(v in c('first_day','no_mpv')){
 z<-fread(file.path(OUT,paste0(v,'_nhanes_scores.csv')));idx<-match(d$SEQN,z$SEQN)
 check(paste(v,'participant linkage'),nrow(z)==44772L&&uniqueN(z$SEQN)==44772L&&!anyNA(idx)&&all(as.character(d$CYCLE)==as.character(z$CYCLE[idx])))
 for(k in c('NM','TB','TC')){
  p<-z[[paste0(k,'_component')]];nobs<-z[[paste0(k,'_input_observed_n')]]
  check(paste(v,k,'observation gate and probability'),identical(is.na(p),nobs==0)&&all(p[!is.na(p)]>=0&p[!is.na(p)]<=1))
 }
 check(paste(v,'Overall gate and range'),identical(is.na(z$Overall_expected_burden),!complete.cases(z[,.(NM_component,TB_component,TC_component)]))&&all(z$Overall_expected_burden[!is.na(z$Overall_expected_burden)]>=0&z$Overall_expected_burden[!is.na(z$Overall_expected_burden)]<=3))
 if(v=='no_mpv')check('TC frozen transfer unchanged',near(z$TC_component[idx][!is.na(d$TC_component)],d$TC_component[!is.na(d$TC_component)]))
}
nf<-fread(file.path(V,'nested_a_primary_final_features.csv'));bi<-fread(file.path(V,'nested_a_primary_boruta_importance.csv'))
check('MPV excluded before selection and final fitting',!'lab_mpv'%in%nf$feature&&!'lab_mpv'%in%bi$feature)
check('No-MPV selected feature counts',identical(as.integer(nf[,.N,by=outcome][order(outcome)]$N),c(11L,7L,22L)))
it<-fread(file.path(V,'nested_a_primary_inner_tuning.csv'));ft<-fread(file.path(V,'nested_a_primary_final_tuning.csv'))
for(o in c('Outcome_NutriMetab','Outcome_TumorBurden')){
 check(paste(o,'nested tuning coverage'),nrow(it[outcome==o])==300L&&nrow(ft[outcome==o])==60L&&uniqueN(it[outcome==o,outer_fold])==5L)
 check(paste(o,'candidate count'),uniqueN(bi[outcome==o,feature])==if(o=='Outcome_NutriMetab')20L else 31L)
}
tcfile<-'catboost_nested_a_primary_Outcome_TreatComp.cbm'
check('TC model file unchanged',unname(tools::md5sum(file.path(V,tcfile)))==unname(tools::md5sum(file.path(REF,'outputs',tcfile))))
oo<-fread(file.path(V,'overall_oof_reproduction_checks.csv'));fr<-fread(file.path(V,'overall_frozen_component_reproduction_checks.csv'))
check('No-MPV outer-fold reconstruction',nrow(oo)==15L&&all(oo$reproduced_within_tolerance)&&max(oo$absolute_auc_difference)<1e-12)
check('No-MPV final component reconstruction',nrow(fr)==3L&&all(fr$within_tolerance)&&max(fr$max_absolute_probability_difference)<1e-12)
for(f in c(file.path(V,'overall_meta_training_oof.csv'),file.path(OUT,'first_day_transfer/overall_meta_training_oof.csv'))){
 a<-fread(f);check(paste(basename(dirname(f)),'OOF participant coverage'),nrow(a)==49082L&&uniqueN(a$patient_key)==49082L&&all(complete.cases(a)))
}
fd<-fread(file.path(OUT,'first_day_transfer/oof_fold_checks.csv'))
check('First-day conditional OOF fold coverage',nrow(fd)==15L&&all(fd$holdout_excluded_from_fit)&&all(fd$train_n+fd$holdout_n==49082L)&&all(fd[,sum(holdout_n),by=component]$V1==49082L))
mp<-fread(file.path(OUT,'missingness/test_predictions.csv'));rp<-fread(file.path(REF,'outputs/nested_a_primary_calendar_test_predictions.csv'))
map<-c(NM='Outcome_NutriMetab',TB='Outcome_TumorBurden',TC='Outcome_TreatComp')
for(k in names(map)){
 a<-mp[component==k];b<-rp[outcome==map[[k]]];idx<-match(a$patient_key,b$patient_key)
 check(paste(k,'primary calendar prediction reproduction'),nrow(a)==11125L&&!anyNA(idx)&&near(a$primary_probability,b$calibrated_probability[idx],1e-12)&&all(a$observed==b$observed[idx]))
 perf<-fread(file.path(OUT,'missingness/hospital_diagnostic_performance.csv'))[component==k & analysis=='Primary model; all inputs observed']
 nmax<-fread(file.path(REF,'outputs/nested_a_primary_final_features.csv'))[outcome==map[[k]],.N]
 check(paste(k,'complete-input subgroup counts'),perf$n==sum(a$observed_inputs==nmax)&&perf$positive==sum(a$observed[a$observed_inputs==nmax]))
}
q<-fread(file.path(OUT,'quartile_counts.csv'));oldq<-unique(fread(file.path(REF,'outputs/supplement_v2/competing_risk_cif_curves.csv'))[,.(score,cause,quartile,n,target_events)])
both<-merge(q,oldq,by=c('score','cause','quartile'),suffixes=c('_new','_old'))
check('Quartile counts reproduce all published CIF strata',nrow(both)==40L&&all(both$n_new==both$n_old)&&all(both$target_events_new==both$target_events_old))
check('Quartile competing event accounting',all(q$n==q$target_events+q$competing_events+q$censored))
qc<-fread(file.path(OUT,'quartile_cutpoints.csv'));qa<-fread(file.path(OUT,'quartile_assignments.csv'))
qd<-d[cipds_complete_mask(d,unname(CIPDS_SCORE_VARS))];raws<-c(NM='NM_component',TB='TB_component',TC='TC_component',OVERALL='Overall_expected_burden',PHENO='PhenoAge_acceleration')
for(k in names(raws)){
 z<-qc[score==k][order(quantile)];ranked<-qd[order(get(raws[[k]]))];cw<-cumsum(ranked$pooled_mec_weight)/sum(ranked$pooled_mec_weight)
 independent<-vapply(c(.25,.5,.75),function(p)ranked[[raws[[k]]]][which(cw>=p)[1]],numeric(1))
 check(paste(k,'weighted empirical quartiles'),near(z$cutpoint_raw,independent,1e-10))
 # Classify with in-memory cutpoints; printed CSV decimals can round a boundary observation.
 a<-qa[score==k];ix<-match(qd$SEQN,a$SEQN);gr<-as.integer(cut(qd[[raws[[k]]]],c(-Inf,independent,Inf),include.lowest=TRUE))
 check(paste(k,'raw-unit group reproduction'),all(gr==a$quartile[ix]))
}
new<-fread(file.path(OUT,'temporal_sensitivity_comparisons.csv'));old<-fread(file.path(REF,'outputs/temporal_combination/temporal_combination_paired_comparisons.csv'))
m<-merge(new[variant=='primary'],old,by='metric',suffixes=c('_new','_old'))
for(k in c('estimate_base','estimate_augmented','paired_difference','ci_lower','ci_upper','p_value','p_holm'))check(paste('Primary temporal reproduction',k),near(m[[paste0(k,'_new')]],m[[paste0(k,'_old')]],1e-8))
for(v in unique(new$variant)){
 a<-new[variant==v]
 check(paste(v,'six paired endpoints and Holm adjustment'),nrow(a)==6L&&near(a$p_holm,p.adjust(a$p_value,'holm')))
 check(paste(v,'paired arithmetic and confidence bounds'),near(a$estimate_augmented-a$estimate_base,a$paired_difference)&&all(a$ci_lower<=a$paired_difference&a$ci_upper>=a$paired_difference)&&near(a$benefit_difference,ifelse(a$higher_is_better,a$paired_difference,-a$paired_difference)))
}
pr<-readRDS(file.path(OUT,'temporal_validation_predictions.rds'));re<-readRDS(file.path(OUT,'temporal_paired_replicates.rds'));fits<-readRDS(file.path(OUT,'temporal_frozen_models.rds'))
check('Bootstrap coverage',nrow(re$replicates)==1000L&&ncol(re$replicates)==60L&&all(is.finite(re$replicates)))
check('Frozen temporal model sample sizes',all(vapply(fits,function(f)f$n,numeric(1))==2579L)&&length(pr$SEQN)==3048L&&sum(pr$event)==1085L)
training_seqn<-d[CYCLE%in%c('2003-2004','2005-2006'),SEQN]
check('Temporal participant separation',length(intersect(training_seqn,pr$SEQN))==0L)
for(k in names(pr$predictions)){
 risk<-pr$predictions[[k]]$risk
 check(paste(k,'risk range and monotonicity'),all(is.finite(risk))&&all(risk>=0&risk<=1)&&all(risk[,-1,drop=FALSE]-risk[,-ncol(risk),drop=FALSE]>=-1e-12))
}
kn<-fread(file.path(OUT,'nonlinear_training_knots.csv'));co<-fread(file.path(OUT,'temporal_phenoage_residualization.csv'));tf<-fread(file.path(OUT,'temporal_training_transforms.csv'))
dom<-d[x$age60 & CYCLE%in%c('2003-2004','2005-2006')];pa<-(dom$PhenoAge-(co$coefficient[co$term=='(Intercept)']+co$coefficient[co$term=='Age']*dom$Age));tr<-tf[variant=='primary'&score=='PHENO'];pz<-(pa-tr$training_mean)/tr$training_sd
for(k in c('Age','PHENO_z_age60'))check(paste(k,'training-only spline knots'),near(kn[variable==k,knot],cipds_weighted_quantile(if(k=='Age')dom$Age else pz,dom$pooled_mec_weight,c(.05,.35,.65,.95))))
check('Spline model basis and component terms',all(paste0('age_ns_',1:3)%in%names(coef(fits$spline_base)))&&all(paste0('pheno_ns_',1:3)%in%names(coef(fits$spline_base)))&&!any(c('Age','PHENO_z_age60')%in%names(coef(fits$spline_base)))&&all(c('NM_z_age60','TB_z_age60','TC_z_age60')%in%names(coef(fits$spline_augmented))))
ass<-fread(file.path(OUT,'variant_mortality_associations.csv'))
check('Association common sample and primary Overall replication',all(ass$n==12094L&ass$events==4787L)&&abs(ass[variant=='primary'&score=='OVERALL',hr]-1.328841)<1e-6)
allchecks<-rbindlist(checks);fwrite(allchecks,file.path(Q,'sensitivity_validation_checks.csv'))
write_json(list(n_checks=nrow(allchecks),n_passed=sum(allchecks$passed),all_passed=all(allchecks$passed),checked_at=as.character(Sys.time())),file.path(Q,'sensitivity_validation_summary.json'),pretty=TRUE,auto_unbox=TRUE)
print(allchecks[passed==FALSE]);cat(sum(allchecks$passed),'of',nrow(allchecks),'checks passed\n');stopifnot(all(allchecks$passed))
