from pathlib import Path
import argparse,os,shutil,subprocess,sys
p=argparse.ArgumentParser(description="Run a purely synthetic demonstration, not manuscript reproduction.")
p.add_argument("--rscript",default=shutil.which("Rscript"))
p.add_argument("--threads",type=int,default=23)
a=p.parse_args()
if not a.rscript: p.error("Supply --rscript; Rscript was not found on PATH.")
root=Path(__file__).resolve().parent
env=os.environ.copy()
for k in list(env):
 if k.startswith("LC_") or k=="LANG":env.pop(k)
env.update(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",R_DATATABLE_NUM_THREADS="1")
print("SYNTHETIC DATA ONLY: this does not reproduce manuscript estimates.",flush=True)
sys.exit(subprocess.call([a.rscript,"--vanilla","demo/run_demo.R",".",str(min(23,max(1,a.threads)))],cwd=root,env=env))

