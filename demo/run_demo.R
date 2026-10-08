# Entirely synthetic demo. No patient records or trained research models are read.
suppressPackageStartupMessages({library(data.table);library(catboost);library(pROC);library(MASS);library(survey);library(survival);library(jsonlite)})
options(survey.lonely.psu="adjust")
a<-commandArgs(TRUE);root<-normalizePath(a[1],winslash="/");out<-file.path(root,"demo_work");dir.create(out,showWarnings=FALSE)
THREADS<-min(23L,max(1L,as.integer(a[2])));setDTthreads(1);SEED<-20261008L;set.seed(SEED)
load_defs<-function(path,wanted){
 found<-character()
 for(e in parse(path,encoding="UTF-8"))if(is.call(e)&&identical(e[[1]],as.name("<-"))&&is.symbol(e[[2]])){
  n<-as.character(e[[2]]);if(n %in% wanted){eval(e,.GlobalEnv);found<-c(found,n)}
 }
 stopifnot(setequal(found,wanted))
}
src<-file.path(root,"analysis/reference/scripts")
load_defs(file.path(src,"07_nested_5x5_boruta_catboost.R"),c("outcome_spec","make_stratified_folds","make_pool","class_weights","fit_catboost","predict_catboost","auc_value","average_precision","boruta_once","tune_inner"))
load_defs(file.path(src,"27_supplement_common.R"),c("cipds_weighted_km_censoring","cipds_weighted_td_auc","cipds_weighted_uno","cipds_weighted_brier","cipds_baseline_hazard","cipds_predict_risk"))
source(file.path(src,"weight_contract.R"),encoding="UTF-8")
scenario<-"SYNTHETIC_DEMO";INNER_K<-2L;OUTER_K<-2L;BORUTA_ITERATIONS<-20L;MODEL_ITERATIONS<-30L
param_grid<-data.table(depth=c(2L,3L),learning_rate=.1,l2_leaf_reg=3)
labs<-c("lab_alb","lab_a_ratio_g","lab_glob","lab_tp","lab_hgb","lab_hct","lab_rbc","lab_rdw","lab_mcv","lab_mch","lab_mchc","lab_wbc","lab_plt","lab_mpv","lab_neu","lab_lym","lab_alp","lab_alt","lab_ast","lab_crea","lab_ggt","lab_tbil","lab_urea","lab_ua")
# Hand-written illustrative centres/scales, not estimated from patient data.
loc<-c(42,1.5,28,70,130,40,4.6,13.5,90,30,335,6,230,10,3.5,1.8,85,25,27,80,35,12,5,300)
spread<-c(4,.2,4,5,15,4,.4,1,5,2,15,1,45,1,1,.5,15,8,8,15,10,3,1,50)
generate<-function(n,prefix){
 latent<-matrix(rnorm(n*3),n,3);d<-data.table(is_synthetic=rep(TRUE,n),synthetic_id=sprintf(paste0(prefix,"%04d"),1:n))
 for(j in seq_along(labs)){
  v<-pmax(loc[j]+spread[j]*(.65*latent[,(j-1)%%3+1]+rnorm(n)),loc[j]*.08)
  v[runif(n)<.07]<-NA_real_;d[[labs[j]]]<-v
 }
 list(data=copy(d),latent=latent)
}
g<-generate(900,"DEMO_H_");h<-g$data;h[,split:=rep(c("train","validation","test"),c(540,180,180))]
outs<-names(outcome_spec);short<-c("NM","TB","TC")
for(k in 1:3)h[[outs[k]]]<-rbinom(nrow(h),1,plogis(-.5+.9*g$latent[,k]))
gs<-generate(1440,"DEMO_S_");s<-gs$data
cycles<-c("1999-2000","2001-2002","2003-2004","2005-2006","2007-2008","2009-2010","2011-2012","2013-2014","2015-2016")
s[,CYCLE:=rep(cycles,each=160)];s[,Age:=runif(.N,45,90)];s[,Sex:=factor(rbinom(.N,1,.5))]
s[,WTMEC2YR:=runif(.N,800,1800)];s[,WTMEC4YR:=runif(.N,900,1700)]
s[,stratum:=paste(CYCLE,rep(1:10,length.out=.N),sep="_")]
s[,psu:=paste(stratum,rep(rep(1:2,each=10),length.out=.N),sep="_")]
s[,stratum:=factor(stratum)];s[,psu:=factor(psu)]
s[,weight:=cipds_pooled_mec_weights(s,9)]
early<-s$CYCLE %in% cycles[1:2]
stopifnot(all.equal(s$weight[early],2*s$WTMEC4YR[early]/9),all.equal(s$weight[!early],s$WTMEC2YR[!early]/9))
observed<-function(v,j){x<-s[[v]];x[is.na(x)]<-loc[j];x}
# Public PhenoAge equation on fictional measurements. Imputed example values are for this demo only.
xb<--19.907-.0336*observed("lab_alb",1)+.0095*observed("lab_crea",20)+.1953*pmax(2,5+rnorm(nrow(s)))+
 .0954*rnorm(nrow(s),-1,1)-.012*rnorm(nrow(s),30,5)+.0268*observed("lab_mcv",9)+
 .3306*observed("lab_rdw",8)+.00188*observed("lab_alp",17)+.0554*observed("lab_wbc",12)+.0804*s$Age
s[,PhenoAge:=141.50+(log(.00553)+log(1.51714)+xb-log(.0076927))/.09165]
et<-rexp(nrow(s),.07*exp(.035*(s$Age-65)+.2*rowSums(gs$latent)));ct<-runif(nrow(s),5,19)
s[,time:=pmin(et,ct)];s[,event:=as.integer(et<=ct)]
tr<-which(h$split=="train");va<-which(h$split=="validation");te<-which(h$split=="test")
oof<-matrix(NA_real_,length(tr),3,dimnames=list(NULL,short));ps<-matrix(NA_real_,nrow(s),3,dimnames=list(NULL,short))
rows<-list();selected<-list();checks<-list()
cal_apply<-function(p,b)plogis(b[1]+b[2]*qlogis(pmin(pmax(p,1e-6),1-1e-6)))
for(k in 1:3){
 outcome<-outs[k];candidate<-setdiff(labs,outcome_spec[[outcome]]$excluded)
 train<-h[tr];folds<-make_stratified_folds(train[[outcome]],2,SEED+k)
 votes<-setNames(integer(length(candidate)),candidate);cat("Synthetic component",short[k],"\n")
 for(f in 1:2){
  i<-which(folds!=f);j<-which(folds==f);stopifnot(!length(intersect(train$synthetic_id[i],train$synthetic_id[j])))
  sel<-boruta_once(train[i,..candidate],train[[outcome]][i],SEED+100*k+f,outcome,f)
  tu<-tune_inner(train[i],sel$selected,outcome,f,SEED+1000*k+f)
  fit<-fit_catboost(train[i,sel$selected,with=FALSE],train[[outcome]][i],as.list(tu$best[,.(depth,learning_rate,l2_leaf_reg)]),SEED+2000*k+f)
  oof[j,k]<-predict_catboost(fit$model,train[j,sel$selected,with=FALSE]);votes[sel$selected]<-votes[sel$selected]+1L
  checks[[length(checks)+1]]<-data.table(is_synthetic=TRUE,component=short[k],fold=f,heldout=length(j),overlap=0)
 }
 features<-names(votes)[votes>=1];stopifnot(!length(intersect(features,outcome_spec[[outcome]]$excluded)))
 tu<-tune_inner(train,features,outcome,0L,SEED+5000*k)
 fit<-fit_catboost(train[,..features],train[[outcome]],as.list(tu$best[,.(depth,learning_rate,l2_leaf_reg)]),SEED+9000*k)
 pv<-predict_catboost(fit$model,h[va,..features]);lp<-qlogis(pmin(pmax(pv,1e-6),1-1e-6))
 cal<-glm(h[[outcome]][va]~lp,family=binomial());b<-unname(coef(cal));stopifnot(all(is.finite(b)))
 oof[,k]<-cal_apply(oof[,k],b);pt<-cal_apply(predict_catboost(fit$model,h[te,..features]),b)
 ps[,k]<-cal_apply(predict_catboost(fit$model,s[,..features]),b)
 ps[rowSums(!is.na(as.matrix(s[,..features])))==0,k]<-NA_real_
 selected[[k]]<-data.table(is_synthetic=TRUE,component=short[k],feature=features)
 rows[[k]]<-data.table(is_synthetic=TRUE,stage="hospital_test_demo",model=short[k],metric=c("AUROC","average_precision","Brier"),
 value=c(auc_value(h[[outcome]][te],pt),average_precision(h[[outcome]][te],pt),mean((h[[outcome]][te]-pt)^2)))
}
stopifnot(all(is.finite(oof)))
tol<-function(p)qlogis(pmin(pmax(p,1e-6),1-1e-6))
z<-apply(oof,2,tol);mu<-colMeans(z);sdv<-apply(z,2,sd)
standard<-function(p)as.data.frame(sweep(sweep(apply(p,2,tol),2,mu,"-"),2,sdv,"/"))
meta<-standard(oof);meta$burden<-ordered(rowSums(as.matrix(h[tr,..outs])),levels=0:3)
stack<-polr(burden~NM+TB+TC,data=meta,method="logistic")
s[,Overall:=as.numeric(predict(stack,standard(ps),type="probs")%*%(0:3))]
for(k in 1:3)s[[short[k]]]<-ps[,k]
stopifnot(all(s$Overall>=0 & s$Overall<=3,na.rm=TRUE))
des<-svydesign(ids=~psu,strata=~stratum,weights=~weight,data=s,nest=TRUE)
fit_hr<-svycoxph(Surv(time,event)~Age+Sex+NM+TB+TC,subset(des,Age>=60))
rows[[4]]<-data.table(is_synthetic=TRUE,stage="survey_association_demo",model="clinical_components",metric=paste0("HR_",names(coef(fit_hr))),value=exp(coef(fit_hr)))
train_idx<-s$Age>=60 & s$CYCLE %in% cycles[3:4];valid_idx<-s$Age>=60 & s$CYCLE %in% cycles[5:6]
freeze<-function(d){
 d<-copy(d[train_idx]);de<-svydesign(ids=~psu,strata=~stratum,weights=~weight,data=d,nest=TRUE)
 co<-coef(svyglm(PhenoAge~Age,design=de));d[,PA:=PhenoAge-(co[1]+co[2]*Age)]
 de<-svydesign(ids=~psu,strata=~stratum,weights=~weight,data=d,nest=TRUE)
 pars<-lapply(c("NM","TB","TC","PA"),function(v)c(mean=as.numeric(coef(svymean(as.formula(paste0("~",v)),de,na.rm=TRUE))),sd=sqrt(as.numeric(svyvar(as.formula(paste0("~",v)),de,na.rm=TRUE)))))
 names(pars)<-c("NM","TB","TC","PA");list(pheno=co,scales=pars)
}
pars<-freeze(s);changed<-copy(s);changed[valid_idx,PhenoAge:=PhenoAge+1000];changed[valid_idx,NM:=.99];changed[valid_idx,event:=1L-event]
stopifnot(identical(pars,freeze(changed)))
for(v in names(pars$scales)){
 raw<-if(v=="PA")s$PhenoAge-(pars$pheno[1]+pars$pheno[2]*s$Age) else s[[v]]
 s[[paste0(v,"_z")]]<-(raw-pars$scales[[v]]["mean"])/pars$scales[[v]]["sd"]
}
train<-s[train_idx];valid<-s[valid_idx]
base<-coxph(Surv(time,event)~Age+Sex+PA_z,data=train,weights=weight/mean(weight),x=TRUE,y=TRUE,model=TRUE)
aug<-coxph(Surv(time,event)~Age+Sex+PA_z+NM_z+TB_z+TC_z,data=train,weights=weight/mean(weight),x=TRUE,y=TRUE,model=TRUE)
for(nm in c("base","aug")){
 pred<-cipds_predict_risk(get(nm),valid,1:10)
 bs<-sapply(1:10,function(t)cipds_weighted_brier(valid$time,valid$event,pred$risk[,t],valid$weight,t))
 rows[[length(rows)+1]]<-data.table(is_synthetic=TRUE,stage="temporal_validation_demo",model=nm,metric=c("Uno_C10","AUC10","IBS_1to10"),
 value=c(cipds_weighted_uno(valid$time,valid$event,pred$lp,valid$weight,10),cipds_weighted_td_auc(valid$time,valid$event,pred$lp,valid$weight,10),sum((head(bs,-1)+tail(bs,-1))/2)/9))
}
p1<-predict(base,valid,type="lp");v2<-copy(valid);v2[,event:=1L-event];stopifnot(identical(p1,predict(base,v2,type="lp")))
metrics<-rbindlist(rows);stopifnot(all(is.finite(metrics$value)))
fwrite(metrics,file.path(out,"demo_metrics.csv"));fwrite(rbindlist(checks),file.path(out,"demo_fold_checks.csv"))
fwrite(rbindlist(selected),file.path(out,"demo_selected_features.csv"))
fwrite(h[1:12,c("is_synthetic","synthetic_id","split",labs[1:6],outs),with=FALSE],file.path(out,"synthetic_preview.csv"))
status<-list(outcome_exclusions=TRUE,disjoint_splits=TRUE,oof_complete=TRUE,validation_calibration=TRUE,
 overall_range=TRUE,corrected_weights=TRUE,training_transforms_invariant_to_validation=TRUE,
 predictions_ignore_validation_labels=TRUE,finite_metrics=TRUE)
write_json(list(scope="SYNTHETIC DEMO ONLY; not manuscript numerical reproduction",seed=SEED,threads=THREADS,
 hospital_n=nrow(h),survey_n=nrow(s),checks=status,all_passed=all(unlist(status))),file.path(out,"demo_validation.json"),pretty=TRUE,auto_unbox=TRUE)
write_json(as.list(installed.packages()[c("data.table","survey","survival","MASS","catboost","pROC","jsonlite"),"Version"]),
 file.path(out,"demo_package_versions.json"),pretty=TRUE,auto_unbox=TRUE)
cat("Synthetic demonstration complete; all checks passed.\n")

