# Controlled reproduction

The public archive deliberately contains no authentic participant data, public-use NHANES records, or trained study models. Reviewers can run the synthetic example without such inputs.

Exact study reproduction requires an institutionally approved external workspace. Depending on the stage, it must contain:

- reference/inputs/nhanes_phase7_batch1.RData: the prepared nine-cycle NHANES object nhanes_valid, with original early-cycle WTMEC4YR and later WTMEC2YR.
- reference/outputs/nhanes_expanded_candidate_matrix.rds and the frozen feature/calibration registry.
- reference/outputs: the privately generated CatBoost component models, ordinal Overall bundle, component and Overall NHANES scores, and prerequisite stage outputs.
- For hospital training: the approved prepared hospital table and candidate registry. These are not recreated from the synthetic preview.
- For sensitivities: the approved first-day and no-MPV component/model artifacts and sensitivity scores.

These filenames document interfaces, not files included in the release. Each source checks its required inputs and study cohort sizes. Do not bypass those gates to make synthetic rows look like real study data. The public source includes the original preparation logic; institution-specific input assembly and full raw-CDC extraction are not certified portable here.

Example plan inspection:

    python run_analysis.py --workspace /secure/project --stage temporal --dry-run

Execution additionally requires --allow-restricted-inputs, an Rscript path if not on PATH, the approved inputs, and prior outputs. The launcher writes only in the separate workspace and refuses to overwrite differing scripts.

Available stages cover hospital A-tier/Overall, primary transport, paired discrimination, apparent increment, temporal comparison, and added NHANES sensitivities. Other source modules are mapped in ANALYSIS_MAP.md and must be executed in dependency order. This is an optional controlled launcher, not an end-to-end reproduction guarantee.

Demo packages tested: R 4.5.1; CatBoost 1.2.5; data.table 1.17.8; pROC 1.18.5; MASS 7.3-65; survey 4.4-8; survival 3.8-3; jsonlite 2.0.0. Full sources additionally use glmnet, mgcv, Hmisc, haven, timeROC and the Python packages in requirements-analysis.txt. See the source manifest for versions of source files. Exact binary reproducibility across operating systems is not claimed.

