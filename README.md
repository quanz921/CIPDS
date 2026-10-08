# CIPDS analysis code and synthetic demonstration

For review of Clinical vulnerability and biological aging in older adults: a cross-cohort study of cancer-trained multidomain laboratory profiles, with S1–S30.

**No real participant records or hospital-trained model objects are included. The synthetic demo does not reproduce the paper's numerical results.**

## Run

Use R 4.5 (or compatible) and Python 3.10+. Install CatBoost R 1.2.5 following official instructions. Run Rscript scripts/install_demo_dependencies.R for other demo packages, then:

    python run_demo.py --rscript /path/to/Rscript --threads 23
    python scripts/check_release.py

On Windows use the full Rscript.exe path. Omit --rscript when it is on PATH. Inspect demo_work/demo_metrics.csv and demo_work/demo_validation.json. All outputs are synthetic. Demo models remain in memory; no model binary is saved.

## Contents

- analysis/reference/scripts: actual analysis source with local paths replaced.
- analysis/sensitivity: missingness, MPV, first-day and nonlinear sensitivity code.
- demo: independent generator and labelled examples.
- docs/ANALYSIS_MAP.md: manuscript mapping and scope.
- docs/SHARING_POLICY.md: privacy boundaries and guidance.
- docs/DATA_AND_CODE_AVAILABILITY.txt: manuscript replacement wording.
- provenance: source hashes and changes, without patient artifacts.

Reviewers can inspect implementation and run a synthetic example. Exact reproduction needs separately authorised hospital inputs/model artifacts and public NHANES preparation. Synthetic results never stand in for real study results.

For controlled work see docs/CONTROLLED_REPRODUCTION.md and run_analysis.py. This is not a data-access mechanism.
The public repository provides source code and synthetic examples only. The copyright holders must select a reuse licence before describing the package as open-source.
