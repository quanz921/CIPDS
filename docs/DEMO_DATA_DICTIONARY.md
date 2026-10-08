# Synthetic example dictionary

All rows have is_synthetic=TRUE and a DEMO_ identifier. No identifier refers to a person.

| Field | Demo meaning |
|---|---|
| synthetic_id | Fictional local row identifier |
| split | Generated train/validation/test assignment |
| lab_* | Random measurements drawn from hand-written constants |
| Outcome_NutriMetab / Outcome_TumorBurden / Outcome_TreatComp | Random binary outcomes generated from hypothetical latent variables |
| CYCLE | Artificial grouping that illustrates the published survey-cycle weighting rule |
| WTMEC2YR / WTMEC4YR | Random positive weights, not NHANES released weights |
| stratum / psu | Artificial survey clusters |
| Age / Sex | Random covariates |
| PhenoAge | Published formula evaluated on fictional measurements |
| time / event | Random follow-up and event indicator |
| NM / TB / TC / Overall | Temporary scores trained exclusively on the generated demo rows |

The demo uses 900 fictional hospital-like rows and 1,440 fictional survey-like rows. Two outer folds, two inner folds, 20 selection trees and 30 model trees are used for speed; these differ from the actual research settings. The clinical demo has age and sex only. Missing inputs for the illustrative PhenoAge equation use hand-written centres; the research sources preserve their actual complete-case rules. No demo performance should be used to interpret the real study.
