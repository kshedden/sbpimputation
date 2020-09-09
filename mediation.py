import sys
sys.path.insert(0, "/afs/umich.edu/user/k/s/kshedden/statsmodels_fork/statsmodels")

import pandas as pd
import statsmodels.api as sm
from statsmodels.imputation.bayes_mi import MI
import os
import patsy
import numpy as np
from data_tools import df
from pca_utils import do_pca
from mimi import mimi
import json
from config import *

n_imp = 20

ref_age = 21

impvar = sys.argv[1]
if impvar not in allowed_controls:
    msg = "Unknown imputation variable: %s" % impvar
    sys.exit(msg)

bp_var = "SBP_MEAN" if sys.argv[2] == "sbp" else "DBP_MEAN"
bp_dir = "sbp" if sys.argv[2] == "sbp" else "dbp"

ndim = int(sys.argv[3])

# Control for father's body size and blood pressure
dadbp = sys.argv[4] == "dadbp"

sex = sys.argv[5]

pat1 = int(sys.argv[6])
pat2 = int(sys.argv[7])

# All output goes here
d = "mixed"
fn = os.path.join(bp_dir, d, "dim_%d" % ndim, "%s_nodadbp_med.txt" % impvar)
if dadbp:
    fn = fn.replace("nodadbp", "dadbp")
if sex != "all":
    fn = fn.replace("med.txt", "%s_med.txt" % sex)
out = open(fn, "w")

m = impvar
if dadbp:
    m += ", with dad BP"
else:
    m += ", with no dad BP"

if sex == "all":
    m += ", all subjects"
elif sex == "female":
    m += ", female subjects only"
else:
    m += ", male subjects only"

out.write(m + "\n")

# load cenv
d1 = "mixed"
d2 = "dim_%d" % ndim
fn = os.path.join(bp_dir, d1, d2, "%s_controlcbs_nodadbp_cen.json" % impvar)
if dadbp:
    fn = fn.replace("nodadbp", "dadbp")
if sex != "all":
    fn = fn.replace("_cen.json", "_%s_cen.json" % sex)
fid = open(fn)
cenv = json.loads(fid.read())
fid.close()

df["Male"] = df["Sex"].replace({"Female": 0, "Male": 1})
df["BMI_cen"] = df.BMI - cenv["BMI"]
df["HT_cen"] = df.HT - cenv["HT"]

df["temp_cen"] = df.Temp - cenv["temp"]
df["age_x"] = (df.Age - 10) / 10

df["ID"] = df["ID"].astype(np.int)

df = df.rename(columns={
    "MomID.Unique": "MomIdUnique",
    "DadID.Unique": "DadIdUnique",
})

varls = [
    "ID", bp_var, "Bamako", "age_x", "Age", "Female",
    "Male", "BMI_cen", "HT_cen", "Year", "School",
    "Wealth_Z_2", "MomIdUnique", "log2T_use_Z", "Breast_Stage_Use_Z",
    "PregMo_Use_cen", "LactMo_Use_cen",
]

varls += ["Mom_BMI", "Mom_Ht_Ave", "Mom_SBP_Z_Res_USE"]

if dadbp:
    varls += ["Dad_BMI", "Dad_Ht_Ave", "Dad_SBP_Z_Res_USE"]

# These are missing for males
ii = (df.Female == 0)
df.loc[ii, "PregMo_Use"] = df.PregMo_Use[ii].fillna(0)
df.loc[ii, "LactMo_Use"] = df.LactMo_Use[ii].fillna(0)

df["LactMo_Use"] = df["LactMo_Use"].fillna(0)

df["LactMo_Use_cen"] = df.LactMo_Use - cenv["LactMo"]
df["PregMo_Use_cen"] = df.PregMo_Use - cenv["PregMo"]

if sex == "female":
	df = df.loc[df.Sex == "Female", :]
elif sex == "male":
	df = df.loc[df.Sex == "Male", :]

log = open(os.path.join(bp_dir, d1, d2, "%s_med.log" % impvar), "w")

mn, cov, proj, vb, vx = do_pca(impvar, ndim, bp_dir)

def do():

    def make_imp(fml):
        imp = MI(
            mimi(impvar, vx, vb, mn, proj, varls, df, log, bp_var, bp_dir),
            sm.OLS,
            formula=fml,
            model_kwds_fn=lambda x: {"data": x},
            burn=0,
            nrep=n_imp,
            skip=0)
        return imp

    rslt, fml = [], []
    for bs_var in "HT_cen", "BMI_cen":
        f = bs_var + " ~ C(Year) + Wealth_Z_2 + Bamako"
        if sex == "all":
            f += " + Female*(age_x + I(age_x**2) + I(age_x**3))"
        else:
            f += " + age_x + I(age_x**2) + I(age_x**3)"
        if sex == "female":
            f += " + Female:Breast_Stage_Use_Z + PregMo_Use_cen + LactMo_Use_cen"
        elif sex == "male":
            f += " + Male:log2T_use_Z "
        f += " + " + " + ".join(vb)
        f += " + Mom_BMI + Mom_Ht_Ave + Mom_SBP_Z_Res_USE"
        if dadbp:
            f += " + Dad_BMI + Dad_Ht_Ave + Dad_SBP_Z_Res_USE"
        fml.append(f)

        imp = make_imp(f)
        r = imp.fit(results_cb=lambda x: x)
        rslt.append(r)

        out.write(r.summary().as_text())
        out.write("\n\n")

    s = np.sqrt(np.diag(cov))
    imp = [make_imp(f).imp for f in fml]
    mnd, mnv, pc = [], [], []
    for k in range(n_imp):

        resid = [None, None]
        scale = [None, None] # unexplained variance
        dz = []
        for j in 0,1:
            imp[j].update()
            dx = imp[j].data
            dx = dx.groupby("ID").head(1)
            dx.loc[:, "age_x"] = (ref_age - 10) / 10
            dz.append(dx)
            xmat = patsy.dmatrix(rslt[j].model.data.design_info, data=dx,
                                 return_type='dataframe')
            resid[j] = dx[bs_var] - np.dot(xmat, rslt[j].params)
            scale[j] = np.mean(resid[j]**2)

        # Get correlated errors for imputed HT and BMI
        r = np.corrcoef(resid[0], resid[1])[0, 1]
        n = dz[0].shape[0]
        ee = []
        for j in 0, 1:
            e1 = np.random.normal(size=n)
            e2 = r*e1 + np.sqrt(1-r**2)*np.random.normal(size=n)
            ee.append([e1, e2])

        yy, pcx = [], []
        patterns = genpat(len(s), s)
        for j in 0, 1:   # HT, BMI
            c = rslt[j].cov_params()
            cr = np.linalg.cholesky(c)
            for jf, f in enumerate([pat1, pat2]):
                pc0 = np.dot(proj.flat, patterns[f][1])
                pcx.append(pc0)
                dz[j].loc[:, impvar + "_pc0"] = pc0
                xmat = patsy.dmatrix(rslt[j].model.data.design_info, data=dz[j],
                                     return_type='dataframe')
                pa = rslt[j].params.copy()
                pa += np.dot(cr, np.random.normal(size=len(pa)))
                y0 = rslt[j].model.predict(params=pa, exog=xmat)
                yy.append(y0)

        pc.append(pcx)
        yd = [yy[1] - yy[0], yy[3] - yy[2]]
        mnd.append([np.mean(yd[0]), np.mean(yd[1])])
        mnv.append([(np.var(yd[0]) + scale[0]) / len(yd[0]),
                    (np.var(yd[1]) + scale[1]) / len(yd[1])])
    mnd = np.asarray(mnd)
    mnv = np.asarray(mnv)
    se = [np.sqrt(np.var(mnd[:, j]) + np.mean(mnv[:, j])) for j in (0,1)]
    mnd = [np.mean(mnd[:, j]) for j in (0,1)]

    return mnd, se, pc[0][1] - pc[0][0], patterns

ym, se, pc, patterns = do()
print(ym)

# Get the model coefficients
d = "mixed"
fn = os.path.join(bp_dir, d, "dim_%d" % ndim, "%s_controlcbs_nodadbp.txt" % impvar)
if sex != "all":
    fn = fn.replace(".txt", "_%s.txt" % sex)
if dadbp:
    fn = fn.replace("nodadbp", "dadbp")
coefs = {}
with open(fn) as inf:
    for line in inf:
        if line.startswith("----------"):
            break
    next(inf)
    next(inf)
    for line in inf:
        if "Var" in line:
            break
        if line.startswith("========="):
            break
        line = line.replace(" ** ", "**")
        v = line.split()
        coefs[v[0]] = (float(v[1]), float(v[2]))


out.write("Childhood %s %s versus %s\n" % (impvar, patterns[pat1][0], patterns[pat2][0]))

cf, cse = [], []
for s in "HT_cen", "BMI_cen", "HT_cen:%s_pc0" % impvar, "BMI_cen:%s_pc0" % impvar:
    cf.append(coefs[s][0])
    cse.append(coefs[s][1])

def bs(x, y, s, t):
    v = [(x+s)*(y+t), (x+s)*(y-t), (x-s)*(y+t), (x-s)*(y-t)]
    return min(v), max(v)

# Indirect effect
z = [0, 0, 0, 0]
z[0] = (cf[0] * ym[0], bs(cf[0], ym[0], cse[0], se[0]))
z[1] = (cf[1] * ym[1], bs(cf[1], ym[1], cse[1], se[1]))
z[2] = (cf[2] * pc, bs(cf[2], pc, cse[2], 0))
z[3] = (cf[3] * pc, bs(cf[3], pc, cse[3], 0))
out.write("Mean change in height at age %d:    %12.3f cm (%.3f)\n" % (ref_age, ym[0], se[0]))
out.write("Mean change in BMI at age %d:       %12.3f kg/m^2 (%.3f)\n" % (ref_age, ym[1], se[1]))

a = z[0][0] + z[2][0]
b = z[0][1][0] + z[2][1][0]
c = z[0][1][1] + z[2][1][1]
out.write("Indirect effect through height:     %12.3f mm/Hg (%.3f,%.3f)\n" % (a, b, c))

a = z[1][0] + z[3][0]
b = z[1][1][0] + z[3][1][0]
c = z[1][1][1] + z[3][1][1]
out.write("Indirect effect through BMI:        %12.3f mm/Hg (%.3f,%.3f)\n" % (a, b, c))

a = sum([x[0] for x in z])
b = sum([x[1][0] for x in z])
c = sum([x[1][1] for x in z])
out.write("Overall indirect effect:            %12.3f mm/Hg (%.3f,%.3f)\n" % (a, b, c))

z = coefs["%s_pc0" % impvar][0] * pc
se = coefs["%s_pc0" % impvar][1] * pc
ci = np.sort([z-se, z+se])
out.write("Direct effect at mean adult size:   %12.3f mm/Hg (%.3f,%.3f)\n" % (z, ci[0], ci[1]))

out.close()
