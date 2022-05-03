"""
Fit models for outcome variables relative to exposures.

Pass the name of the outcome variable

> python sbp.py impvar [controlcbs|nocontrolcbs]
"""

import sys
sys.path.insert(
    0, "/afs/umich.edu/user/k/s/kshedden/statsmodels_fork/statsmodels")

import pandas as pd
import statsmodels.api as sm
from statsmodels.imputation.bayes_mi import MI
import os
import numpy as np
from data_tools import df
from pca_utils import do_pca
from mimi import mimi
import json
from config import *

impvar = sys.argv[1]
if impvar not in allowed_controls:
    msg = "Unknown imputation variable: %s" % impvar
    sys.exit(msg)

bp_var = "SBP_MEAN" if sys.argv[2] == "sbp" else "DBP_MEAN"
bp_dir = "sbp" if sys.argv[2] == "sbp" else "dbp"

ndim = int(sys.argv[3])

# Control for current body size
controlcbs = sys.argv[4] == "controlcbs"

# Control for father's body size and blood pressure
dadbp = sys.argv[5] == "dadbp"

# Save means before centering
cenv = {}

df["Male"] = df["Sex"].replace({"Female": 0, "Male": 1})
cenv["BMI"] = df.BMI.mean()
df["BMI_cen"] = df.BMI - cenv["BMI"]
cenv["HT"] = df.HT.mean()
df["HT_cen"] = df.HT - cenv["HT"]

df.HT.mean()
cenv["temp"] = df.Temp.mean()
df["temp_cen"] = df.Temp - cenv["temp"]
cenv["age"] = df.Age.mean()
df["age_cen"] = df.Age - cenv["age"]
df["age_x"] = (df.Age - 10) / 10

df["ID"] = df["ID"].astype(np.int)

df.loc[df.Bamako==1, "Village"] = 50
df = df.loc[pd.notnull(df["Village"]), :]
df["Village"] = df["Village"].astype(np.int)

df = df.rename(columns={
    "MomID.Unique": "MomIdUnique",
    "DadID.Unique": "DadIdUnique",
})

varls = [
    "ID", bp_var, "age_cen", "age_x", "Age", "Female",
    "Male", "BMI_cen", "HT_cen", "Year", "lognummeas", "temp_cen", "School",
    "Wealth_Z_2", "MomIdUnique", "log2T_use_Z", "Breast_Stage_Use_Z",
    "PregMo_Use_cen", "lactMo_Use_cen", "Village",
]

varls += ["Mom_BMI", "Mom_Ht_Ave", "Mom_SBP_Z_Res_USE"]

if dadbp:
    varls += ["Dad_BMI", "Dad_Ht_Ave", "Dad_SBP_Z_Res_USE"]

# These are missing for males
ii = (df.Female == 0)
df.loc[ii, "PregMo_Use"] = df.PregMo_Use[ii].fillna(0)
df.loc[ii, "lactMo_Use"] = df.lactMo_Use[ii].fillna(0)

df["lactMo_Use"] = df["lactMo_Use"].fillna(0)

cenv["lactMo"] = df.lactMo_Use.mean()
df["lactMo_Use_cen"] = df.lactMo_Use - cenv["lactMo"]
cenv["PregMo"] = df.PregMo_Use.mean()
df["PregMo_Use_cen"] = df.PregMo_Use - cenv["PregMo"]

# Save the means to a file
d1 = "mixed"
d2 = "dim_%d" % ndim
fn = os.path.join(bp_dir, d1, d2, "%s_nocontrolcbs_nodadbp_cen.json" % impvar)
if controlcbs:
    fn = fn.replace("nocontrolcbs", "controlcbs")
if dadbp:
    fn = fn.replace("nodadbp", "dadbp")
fid = open(fn, "w")
fid.write(json.dumps(cenv))
fid.close()

log = open(os.path.join(bp_dir, d1, d2, "%s.log" % impvar), "w")

mn, cov, proj, vb, vx = do_pca(impvar, ndim, bp_dir)

fml = bp_var + " ~ Female*(age_x + I(age_x**2) + I(age_x**3) + lognummeas) + C(Year) + temp_cen + School + Wealth_Z_2 + Female:Breast_Stage_Use_Z + Male:log2T_use_Z + PregMo_Use_cen + I(PregMo_Use_cen**2) + lactMo_Use_cen"

fml += " + C(Village, Treatment(reference=50)) + Mom_BMI + Mom_Ht_Ave + Mom_SBP_Z_Res_USE"

if dadbp:
    fml += " + Dad_BMI + Dad_Ht_Ave + Dad_SBP_Z_Res_USE"

# Controls for body dimensions at time of SBP
if controlcbs:
    fml += " + (HT_cen + BMI_cen)*" + "(" + " + ".join(vb) + ")"
else:
    fml += " + " + " + ".join(vb)

d = "mixed"
fn = os.path.join(bp_dir, d, "dim_%d" % ndim, "%s_nocontrolcbs_nodadbp.txt" % impvar)
if controlcbs:
    fn = fn.replace("nocontrolcbs", "controlcbs")
if dadbp:
    fn = fn.replace("nodadbp", "dadbp")
out = open(fn, "w")

vcf = {"Id": "0 + C(ID)", "Id*age": "0 + C(ID)*age_cen"}

#vcf = {"Id": "0 + C(ID)"}

def model_kwds_fn(x):
    return {"groups": "MomIdUnique", "data": x, "re_formula": "1", "vc_formula": vcf}

def fit_kwds_fn(x):
    return {"method": "lbfgs", "reml": False}


imp = MI(
    mimi(impvar, vx, vb, mn, proj, varls, df, log, bp_var, bp_dir),
    sm.MixedLM,
    model_args_fn=None,
    formula=fml,
    model_kwds_fn=model_kwds_fn,
    fit_kwds=fit_kwds_fn,
    burn=0,
    nrep=20,
    skip=0)

if (impvar == "BMI") and (ndim == 1) and controlcbs:

    import json
    f = open("centering.json")
    centering = json.load(f)
    f.close()

    # Create a table 1, based on counting people
    m = mimi(impvar, vx, vb, mn, proj, varls, df, log, bp_var, bp_dir)
    m.update()
    x = m.data.copy()
    x["temp"] = x["temp_cen"] + cenv["temp"]
    x["BMI"] = x["BMI_cen"] + cenv["BMI"]
    x["HT"] = x["HT_cen"] + cenv["HT"]
    x["NumMeas"] = 10**(x.lognummeas)
    ms = centering["log2T_use"]
    x["Testosterone"] = 2**(x.log2T_use_Z * ms[1] + ms[0])
    ms = centering["Breast_Stage_Use"]
    x["Breast_Stage"] = x.Breast_Stage_Use_Z * ms[1] + ms[0]
    x = x.sort_values(by=["ID", "Age"])
    first = x.groupby("ID").head(1)
    last = x.groupby("ID").tail(1)
    stats = ["Age", "SBP_MEAN", "HT", "BMI", "School", "NumMeas", "temp",
             "Wealth_Z_2", "Breast_Stage", "Testosterone"]
    def q10(x):
        return x.quantile(0.1)
    def q90(x):
        return x.quantile(0.9)
    stats = {v: [len, np.mean, q10, q90] for v in stats}
    astats = []
    for dz in first, last:
        for female in 0, 1:
            dx = dz.loc[dz.Female == female, :]
            a = dx.agg(stats)
            a = a.T
            a["Female"] = female
            a["Visit"] = "first" if dz is first else "last"
            astats.append(a)
    astats = pd.concat(astats, axis=0)
    astats = astats.rename(columns={"mean": "Mean", "len": "N"})
    sname = "%s_table1.csv" % bp_var
    sname = sname.lower()
    sname = sname.replace("mean_", "")
    astats.to_csv(sname)

    # Table based on observations
    y = x.groupby("Female").ID.agg(len)

rslt = imp.fit(results_cb=lambda x: x)

mm = rslt.results[0].model
nobs = sum([x.shape[1] for x in mm.exog_vc.mats[0]])

ic = [x.llf for x in rslt.results]

sca = [x.scale for x in rslt.results]

out.write("%s\n" % impvar)
out.write("%d distinct subjects\n" % nobs)
out.write("%d distinct mothers\n" % mm.n_groups)
out.write("mean IC %f\n" % np.mean(ic))
out.write("mean scale %f\n" % np.mean(sca))
out.write(rslt.summary().as_text())

# Save the coefficient trajectories
if controlcbs:
    vb1 = ["HT_cen", "BMI_cen"] + vb + ["HT_cen:%s" % y for y in vb] + ["BMI_cen:%s" % y for y in vb]
else:
    vb1 = vb
mp = rslt.params
cm = rslt.cov_params()
xn = rslt.model.exog_names
mp = pd.Series(mp[0:len(xn)], xn)
mp = mp[vb1]
cm = pd.DataFrame(cm[0:len(xn), 0:len(xn)], index=xn, columns=xn)
cm = cm.loc[vb1, vb1]

out.write("BEGIN-PROJECTION\n")
for k in range(proj.shape[0]):
    out.write(",".join([str(x) for x in proj[k, :]]) + "\n")
out.write("--\n")

# Save the average within-imputation variance of the growth variable at each age.
for v in np.diag(cov):
    out.write("%f\n" % v)
out.write("--\n")

# Save some of the estimated model parameters that we need later
out.write(mp.to_csv())
out.write("--\n")

# Save the covariance (sampling uncertainty) for the key model parameters
out.write(cm.to_csv())
out.write("END-PROJECTION\n\n")

out.close()
log.close()
