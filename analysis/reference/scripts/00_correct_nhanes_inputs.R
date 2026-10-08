suppressPackageStartupMessages({library(data.table); library(jsonlite)})
root <- Sys.getenv("CIPDS_PACKAGE_DIR", "cipds_corrected_20260908")
env <- new.env()
load("cipds_rebuild_20260831/inputs/nhanes_phase7_batch1.RData", envir=env)
d <- as.data.table(env$nhanes_valid)
stopifnot(nrow(d)==44772L, uniqueN(d$SEQN)==44772L)
original <- copy(d)
coalesce_num <- function(a,b) {z<-as.numeric(a); z[is.na(z)]<-as.numeric(b)[is.na(z)];z}
d[, Albumin := coalesce_num(LBDSALSI, LBXSAL*10)]
d[, Creatinine_uncalibrated_mgdL := coalesce_num(LBXSCR, LBDSCRSI/88.4)]
d[, Creatinine_calibrated_mgdL := Creatinine_uncalibrated_mgdL]
# CDC LAB18 and BIOPRO_D: calibrate once from the archived original values.
d[CYCLE=="1999-2000", Creatinine_calibrated_mgdL := 0.147+1.013*Creatinine_uncalibrated_mgdL]
d[CYCLE=="2005-2006", Creatinine_calibrated_mgdL := -0.016+0.978*Creatinine_uncalibrated_mgdL]
d[, Creatinine := Creatinine_calibrated_mgdL*88.4]
d[, Alkaline_Phosphatase := coalesce_num(LBXSAPSI, LBDSAPSI)]
d[, Glucose := coalesce_num(LBDSGLSI, LBXSGL/18.01559)]
d[, PA_Albumin_gdL := Albumin/10]
d[, PA_Albumin_gL := Albumin]
d[, PA_Creatinine_mgdL := Creatinine_calibrated_mgdL]
d[, PA_Creatinine_umolL := Creatinine]
d[, PA_Glucose_mgdL := coalesce_num(LBXSGL, Glucose*18.01559)]
d[, PA_Glucose_mmolL := Glucose]
d[, PA_ALP_UL := Alkaline_Phosphatase]
d[, PA_CRP_mgdL := coalesce_num(LBXCRP, LBXHSCRP/10)]
stopifnot(!any(d$PA_CRP_mgdL<=0,na.rm=TRUE))
d[, PA_xb := -19.907 -0.0336*PA_Albumin_gL +0.0095*PA_Creatinine_umolL +
  0.1953*PA_Glucose_mmolL +0.0954*log(PA_CRP_mgdL) -0.012*PA_Lymph_pct +
  0.0268*PA_MCV +0.3306*PA_RDW +0.00188*PA_ALP_UL +0.0554*PA_WBC_1000+0.0804*Age]
# Published 2019 corrected equation, algebraically evaluated in log space.
d[, PhenoAge := 141.50+(log(0.00553)+log(1.51714)+PA_xb-log(0.0076927))/0.09165]
d[, PA_mort_prob := -expm1(-exp(PA_xb)*1.51714/0.0076927)]
# The analysis scripts calculate survey-weighted regression residuals within each age domain.
d[, PhenoAge_Accel := PhenoAge - predict(lm(PhenoAge ~ Age, data=d, weights=WTMEC2YR), newdata=d)]
d[, Follow_Up_Months := fifelse(PERMTH_EXM==0,0.5,as.numeric(PERMTH_EXM))]
d[, Follow_Up_Years := Follow_Up_Months/12]
d[, Diabetes := fifelse(DIQ010==1,1L,fifelse(DIQ010 %in% c(2,3),0L,NA_integer_))]
d[, Hypertension := fifelse(BPQ020==1,1L,fifelse(BPQ020==2,0L,NA_integer_))]
cv <- as.matrix(d[,.(MCQ160B,MCQ160C,MCQ160D,MCQ160E,MCQ160F)])
yes <- rowSums(cv==1,na.rm=TRUE)>0
no <- rowSums(cv==2,na.rm=TRUE)==ncol(cv)
d[, CVD := fifelse(yes,1L,fifelse(no,0L,NA_integer_))]
d[, CVD_mortality_eligible := CYCLE!="2015-2016"]
# Preserve the historical ALQ101 threshold definition; explain it explicitly in tables.
nhanes_valid <- as.data.frame(d)
save(nhanes_valid,file=file.path(root,"inputs/nhanes_phase7_batch1.RData"),compress="gzip")
mapped <- as.data.table(readRDS(file.path(root,"outputs/nhanes_expanded_candidate_matrix.rds")))
mi <- match(mapped$SEQN,d$SEQN); stopifnot(!anyNA(mi))
mapped[, lab_crea := d$Creatinine[mi]]
mapped[, lab_glu := d$Glucose[mi]]
saveRDS(mapped,file.path(root,"outputs/nhanes_expanded_candidate_matrix.rds"),compress="gzip")
vars <- c("Age","Albumin","Creatinine","Glucose","Alkaline_Phosphatase","PA_CRP_mgdL","PhenoAge","Follow_Up_Years","Diabetes","Hypertension","CVD")
changes <- rbindlist(lapply(vars,function(v) data.table(variable=v,
  old_nonmissing=sum(!is.na(original[[v]])),new_nonmissing=sum(!is.na(d[[v]])),
  changed_value=sum(abs(original[[v]]-d[[v]])>1e-10,na.rm=TRUE),
  missingness_changed=sum(is.na(original[[v]])!=is.na(d[[v]])))))
fwrite(changes,file.path(root,"qa/input_changes.csv"))
cycle <- d[,.(n=.N,age60=sum(Age>=60),pheno_complete=sum(is.finite(PhenoAge)),
  pheno_complete_age60=sum(is.finite(PhenoAge)&Age>=60),exam_zero_months=sum(PERMTH_EXM==0)),by=CYCLE]
fwrite(cycle,file.path(root,"qa/corrected_input_by_cycle.csv"))
verify <- d[,.(SEQN,CYCLE,Age,PA_Albumin_gL,PA_Creatinine_umolL,PA_Glucose_mmolL,
 PA_CRP_mgdL,PA_Lymph_pct,PA_MCV,PA_RDW,PA_ALP_UL,PA_WBC_1000,PA_xb,PhenoAge,
 PERMTH_INT,PERMTH_EXM,Follow_Up_Years,DIQ010,Diabetes,BPQ020,Hypertension,
 MCQ160B,MCQ160C,MCQ160D,MCQ160E,MCQ160F,CVD,Creatinine_uncalibrated_mgdL,Creatinine_calibrated_mgdL)]
fwrite(verify,file.path(root,"qa/corrected_input_verification.csv"))
write_json(list(phenoage_equation="Liu et al 2019 correction; stable algebraic form", glucose="linear mmol/L",creatinine="CDC calibrated umol/L",albumin="g/L",crp="natural logarithm of mg/dL; no arbitrary lower floor",followup="PERMTH_EXM; recorded zero months assigned 0.5 month",cvd_mortality="heart disease and cerebrovascular death in 1999-2014",unknown_covariates="7/9 and incomplete negative composite responses are missing",alcohol="ever >=12 drinks in any one year (ALQ101)",source_n=nrow(d)),file.path(root,"qa/input_definition_manifest.json"),auto_unbox=TRUE,pretty=TRUE)
print(changes);print(cycle)
