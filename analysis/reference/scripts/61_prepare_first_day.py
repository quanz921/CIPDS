from pathlib import Path
import pandas as pd
import numpy as np
import json
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs'/'first_day_sensitivity'
OUT.mkdir(exist_ok=True)
release=Path('hospital_source')
features=pd.read_csv(ROOT/'outputs/nested_a_primary_final_features.csv').feature.unique().tolist()
core=pd.read_csv(ROOT/'outputs/hospital_patient_model_ready_129.csv',low_memory=False)
entry=pd.to_datetime(core.first_observed_test_date,errors='raise').dt.normalize()
dates=pd.DataFrame({'patient_key':core.patient_key,'entry_date':entry,'split':core.split_calendar_entry})
lab=pd.read_parquet(release/'162_patient_lab_summary_long_v1.parquet',
    columns=['patient_key','final_variable_id','first_numeric_value','first_numeric_datetime'],
    filters=[('final_variable_id','in',features)])
lab=lab[lab.first_numeric_value.notna()].merge(dates,on='patient_key',validate='many_to_one')
lab['date']=pd.to_datetime(lab.first_numeric_datetime,errors='raise').dt.normalize()
lab['lag_days']=(lab.date-lab.entry_date).dt.days
assert lab.lag_days.notna().all() and (lab.lag_days>=0).all()
original=lab.pivot(index='patient_key',columns='final_variable_id',values='first_numeric_value').reindex(core.patient_key)
for f in features:
    assert np.allclose(core[f].to_numpy(),original[f].to_numpy(),equal_nan=True,atol=1e-7,rtol=0)
same=lab[lab.lag_days==0].pivot(index='patient_key',columns='final_variable_id',values='first_numeric_value').reindex(core.patient_key)
keep=['patient_key','first_observed_test_date','split_calendar_entry','Outcome_NutriMetab','Outcome_TumorBurden','Outcome_TreatComp']
result=core[keep].copy()
for f in features:result[f]=same[f].to_numpy()
result.to_csv(OUT/'hospital_first_day_inputs.csv',index=False)
rows=[]
ends={'train':'2024-07-18','validation':'2025-02-03','test':'2025-07-22'}
for split,end in ends.items():
    sub=lab[lab.split==split]
    pp=sub.groupby('patient_key').agg(first_input=('date','min'),last_input=('date','max'),max_lag=('lag_days','max'))
    span=(pp.last_input-pp.first_input).dt.days
    rows.append({'split':split,'patients':int((core.split_calendar_entry==split).sum()),
        'patients_with_selected_input':len(pp),'any_later_than_entry':int((pp.max_lag>0).sum()),
        'any_after_split_end':int((pp.last_input>pd.Timestamp(end)).sum()),
        'input_span_median':float(span.median()),'input_span_q25':float(span.quantile(.25)),
        'input_span_q75':float(span.quantile(.75)),'input_span_p90':float(span.quantile(.9)),
        'input_span_max':int(span.max()),'first_day_values':int((sub.lag_days==0).sum()),
        'later_values_masked':int((sub.lag_days>0).sum()),'first_day_after_split_end':0})
pd.DataFrame(rows).to_csv(OUT/'sampling_window_audit.csv',index=False)
assert not any((lab.loc[lab.lag_days==0,'date']>lab.loc[lab.lag_days==0,'split'].map(ends).pipe(pd.to_datetime)))
(OUT/'preparation_qa.json').write_text(json.dumps({'n':len(result),'selected_features':len(features),
    'all_original_nonmissing_values_verified':len(lab),'sampling_dates_missing':0,
    'sensitivity_window':'first recorded laboratory test calendar day',
    'all_sensitivity_inputs_within_entry_day_and_split':True,
    'feature_sets_and_hyperparameters':'fixed from primary analysis; no new feature search',
    'labels':'unchanged patient-ever documented flags'},indent=2),encoding='utf-8')
print(json.dumps(rows,indent=2))
