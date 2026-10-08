suppressPackageStartupMessages({library(data.table);library(catboost);library(pROC);library(jsonlite)})
setDTthreads(1)
root<-Sys.getenv('CIPDS_PACKAGE_DIR','cipds_revision_20260915')
out<-file.path(root,'outputs/first_day_sensitivity')
d<-fread(file.path(out,'hospital_first_day_inputs.csv'))
original<-fread(file.path(root,'outputs/hospital_patient_model_ready_129.csv'))
features<-fread(file.path(root,'outputs/nested_a_primary_final_features.csv'))
registry<-fread(file.path(root,'outputs/nested_a_primary_calibration_and_model_registry.csv'))
stopifnot(identical(d$patient_key,original$patient_key))
outcomes<-c('Outcome_NutriMetab','Outcome_TumorBurden','Outcome_TreatComp')
make_pool<-function(z,selected,y=NULL) {
 x<-as.matrix(z[,..selected]);storage.mode(x)<-'double'
 if(is.null(y))catboost.load_pool(x,feature_names=as.list(selected)) else
   catboost.load_pool(x,label=as.numeric(y),feature_names=as.list(selected))
}
ap<-function(y,p) {
 o<-order(p,decreasing=TRUE); yy<-y[o]; pp<-p[o]
 ends<-c(which(diff(pp)!=0),length(pp));tp<-cumsum(yy)[ends]
 sum(diff(c(0,tp))/sum(y)*(tp/ends))
}
eps<-1e-6
rows<-list();calrows<-list();prediction_rows<-list();checks<-list()
for(oi in seq_along(outcomes)) {
 outcome<-outcomes[oi]; code_name<-outcome
 selected<-features[outcome==code_name,feature]
 reg<-registry[outcome==code_name]
 stopifnot(!anyDuplicated(selected),length(selected)==c(12L,24L,7L)[oi],nrow(reg)==1L)
 train<-d[split_calendar_entry=='train' & !is.na(get(outcome))]
 val<-d[split_calendar_entry=='validation' & !is.na(get(outcome))]
 test<-d[split_calendar_entry=='test' & !is.na(get(outcome))]
 y<-train[[outcome]]
 params<-list(loss_function='Logloss',eval_metric='AUC',iterations=as.integer(reg$iterations),
  depth=as.integer(reg$depth),learning_rate=as.numeric(reg$learning_rate),
  l2_leaf_reg=as.numeric(reg$l2_leaf_reg),random_seed=20260831L+oi*100000L+999L,
  thread_count=7L,class_weights=c(1,sum(y==0)/sum(y==1)),logging_level='Silent',allow_writing_files=FALSE)
 fit<-catboost.train(make_pool(train,selected,y),params=params)
 catboost.save_model(fit,file.path(out,paste0(outcome,'_first_day.cbm')))
 pv<-as.numeric(catboost.predict(fit,make_pool(val,selected),prediction_type='Probability'))
 val_lp<-qlogis(pmin(pmax(pv,eps),1-eps));vy<-val[[outcome]]
 cal<-glm(vy~val_lp,family=binomial())
 coefcal<-unname(coef(cal))
 pt<-as.numeric(catboost.predict(fit,make_pool(test,selected),prediction_type='Probability'))
 ps<-plogis(coefcal[1]+coefcal[2]*qlogis(pmin(pmax(pt,eps),1-eps)))
 original_test<-original[match(test$patient_key,original$patient_key)]
 oldfit<-catboost.load_model(file.path(root,'outputs',reg$model_file))
 oldraw<-as.numeric(catboost.predict(oldfit,make_pool(original_test,selected),prediction_type='Probability'))
 po<-plogis(reg$intercept+reg$slope*qlogis(pmin(pmax(oldraw,eps),1-eps)))
 yy<-test[[outcome]]
 ro<-roc(yy,po,direction='<',quiet=TRUE);rs<-roc(yy,ps,direction='<',quiet=TRUE)
 ao<-as.numeric(auc(ro));as<-as.numeric(auc(rs));ci<-as.numeric(ci.auc(rs,method='delong'))
 diff<-as-ao;se<-sqrt(max(0,as.numeric(var(rs)+var(ro)-2*cov(rs,ro))))
 set.seed(2026091500L+oi)
 ip<-which(yy==1);ineg<-which(yy==0)
 apboot<-replicate(1000,{idx<-c(sample(ip,length(ip),TRUE),sample(ineg,length(ineg),TRUE));
   c(original=ap(yy[idx],po[idx]),first_day=ap(yy[idx],ps[idx]))})
 apci<-quantile(apboot['first_day',],c(.025,.975))
 lp<-qlogis(pmin(pmax(ps,eps),1-eps));testcal<-glm(yy~lp,family=binomial())
 rows[[oi]]<-data.table(outcome,train_n=nrow(train),validation_n=nrow(val),test_n=nrow(test),
  test_positive=sum(yy),original_auc=ao,first_day_auc=as,first_day_auc_lower=ci[1],first_day_auc_upper=ci[3],
  paired_auc_difference=diff,paired_auc_lower=diff-qnorm(.975)*se,paired_auc_upper=diff+qnorm(.975)*se,
  paired_auc_p=2*pnorm(-abs(diff/se)),original_ap=ap(yy,po),first_day_ap=ap(yy,ps),
  first_day_ap_lower=apci[1],first_day_ap_upper=apci[2],
  original_brier=mean((yy-po)^2),first_day_brier=mean((yy-ps)^2),
  first_day_test_calibration_intercept=coef(testcal)[1],first_day_test_calibration_slope=coef(testcal)[2],
  first_day_test_all_inputs_missing=sum(rowSums(!is.na(test[,..selected]))==0))
 calrows[[oi]]<-data.table(outcome,intercept=coefcal[1],slope=coefcal[2],depth=reg$depth,
  learning_rate=reg$learning_rate,l2_leaf_reg=reg$l2_leaf_reg,iterations=reg$iterations,features=length(selected))
 prediction_rows[[oi]]<-data.table(patient_key=test$patient_key,outcome,observed=yy,original_probability=po,first_day_probability=ps)
 cat(outcome,' completed; AUC=',as,' original=',ao,'\n');flush.console()
}
ans<-rbindlist(rows);ans[,paired_auc_p_holm:=p.adjust(paired_auc_p,'holm')]
fwrite(ans,file.path(out,'first_day_calendar_test_comparison.csv'),bom=TRUE)
fwrite(rbindlist(calrows),file.path(out,'first_day_model_registry.csv'))
fwrite(rbindlist(prediction_rows),file.path(out,'first_day_test_predictions.csv'))
write_json(list(models=3,feature_selection='Primary feature sets retained',
 hyperparameters='Primary final parameters and seeds retained',
 training='First-day inputs, original training patients, outcome-specific eligibility',
 calibration='Refitted only in original calendar-validation cohort',
 testing='Original independent patient calendar-test cohort; first-day inputs',
 auc_interval='DeLong',auc_difference='Paired DeLong; Holm across three outcomes',
 ap_interval='1000 outcome-stratified patient bootstrap percentile intervals',
 hospital_primary_models_modified=FALSE,nhanes_primary_models_modified=FALSE),
 file.path(out,'sensitivity_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
print(ans)
