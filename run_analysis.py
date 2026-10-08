"""Optional launcher for an independently authorised external workspace."""
from pathlib import Path
import argparse,json,os,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parent
STAGES={
 "hospital":[("reference","07_nested_5x5_boruta_catboost.R"),("reference","14_train_overall_stacked_model.R")],
 "primary":[("reference","10_score_nhanes_nested_survey.R"),("reference","12_nhanes_older_adult_primary.R"),("reference","15_nhanes_overall_model_primary.R")],
 "paired":[("reference","19_nhanes_paired_discrimination.R")],
 "incremental":[("reference","30_incremental_prediction.R")],
 "temporal":[("reference","62_temporal_combination_increment.R")],
 "sensitivities":[("sensitivity","91_nhanes_sensitivities.R")],
}
p=argparse.ArgumentParser(description=__doc__)
p.add_argument("--workspace",required=True,type=Path)
p.add_argument("--stage",choices=STAGES,required=True)
p.add_argument("--rscript",default=shutil.which("Rscript"))
p.add_argument("--threads",default=23,type=int)
p.add_argument("--dry-run",action="store_true")
p.add_argument("--allow-restricted-inputs",action="store_true")
a=p.parse_args();work=a.workspace.resolve()
if work==ROOT or ROOT in work.parents or work in ROOT.parents:
 p.error("The controlled workspace must be separate from the public repository.")
plan={"stage":a.stage,"scripts":STAGES[a.stage],"threads":min(23,max(1,a.threads)),
 "data_in_archive":False,"full_reproduction_verified":False,
 "notice":"Prepared authorised inputs and prerequisite outputs must exist; this does not grant data access."}
if a.dry_run:print(json.dumps(plan,indent=2));sys.exit(0)
if not a.allow_restricted_inputs:p.error("Explicit --allow-restricted-inputs is required for an authorised controlled workspace.")
if not a.rscript:p.error("Rscript is required.")
if not work.is_dir():p.error("Provide an existing authorised workspace; this launcher does not create data.")
for unit in ("reference","sensitivity"):
 for source in (ROOT/"analysis"/unit).rglob("*"):
  if not source.is_file():continue
  target=work/unit/source.relative_to(ROOT/"analysis"/unit)
  if target.exists() and target.read_bytes()!=source.read_bytes():
   p.error("Different code already exists at "+str(target)+"; reconcile versions explicitly.")
  target.parent.mkdir(parents=True,exist_ok=True)
  if not target.exists():shutil.copyfile(source,target)
env=os.environ.copy()
for k in list(env):
 if k.startswith("LC_") or k=="LANG":env.pop(k)
env.update(CIPDS_REFERENCE_DIR="reference",CIPDS_SENSITIVITY_DIR="sensitivity",
 CIPDS_THREADS=str(min(23,max(1,a.threads))),CIPDS_SUPPLEMENT_REPLICATES="1000",
 CIPDS_DISCRIMINATION_REPLICATES="1000",OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
for unit,name in STAGES[a.stage]:
 env["CIPDS_PACKAGE_DIR"]="reference" if unit=="sensitivity" else unit
 result=subprocess.run([a.rscript,"--vanilla",f"{unit}/scripts/{name}"],cwd=work,env=env)
 if result.returncode:sys.exit(result.returncode)

