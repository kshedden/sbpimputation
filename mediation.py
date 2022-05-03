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

np.random.seed(1234)

n_imp = 20

ref_age = 21

if len(sys.argv) != 8:
    1/0

impvar = sys.argv[1]
if impvar not in allowed_controls:
    msg = "Unknown imputation variable: %s" % impvar
    sys.exit(msg)

bp_var = "SBP_MEAN" if sys.argv[2] == "sbp" else "DBP_MEAN"
bp_dir = "sbp" if sys.argv[2] == "sbp" else "dbp"

ndim = int(sys.argv[3])

sex = sys.argv[4]

pat1 = int(sys.argv[5])
pat2 = int(sys.argv[6])

adj = int(sys.argv[7])

# All output goes here
d = "mixed"
fn = os.path.join(bp_dir, d, "dim_%d" % ndim, "%s_adj%d_%s_med.txt" % (impvar, adj, sex))
out = open(fn, "w")

msg = impvar

if sex == "female":
    msg += ", female subjects only"
elif sex == "male":
    msg += ", male subjects only"
else:
    1/0

msg += ", adjustment=%d" % adj

msg = "==== " + msg + " ===="
b = (76 - len(msg)) // 2
msg = " "*b + msg
out.write(msg + "\n\n")

# load cenv
d1 = "mixed"
d2 = "dim_%d" % ndim
fn = os.path.join(bp_dir, d1, d2, "%s_adj%d_%s_cen.json" % (impvar, adj, sex))
fid = open(fn)
cenv = json.loads(fid.read())
fid.close()

df["Male"] = df["Sex"].replace({"Female": 0, "Male": 1})
df["BMI_cen"] = df.BMI - cenv["BMI"]
df["HT_cen"] = df.HT - cenv["HT"]

df["temp_cen"] = df.Temp - cenv["temp"]
df["age_x"] = (df.Age - 10) / 10

df["PregMo_Use"] = df["PregMo_Use"].fillna(0)
ii = pd.isnull(df.no_lactmo_info) & pd.isnull(df.lactmo_Use)
df.loc[ii, "lactmo_Use"] = 0
ii = (df.Female == 0)
df.loc[ii, "lactmo_Use"] = 0

df["ID"] = df["ID"].astype(np.int)

if sex == "female":
    df["PregMo_Use_cen"] = df["PregMo_Use"] - cenv["PregMo"]
    df["lactmo_Use_cen"] = df["lactmo_Use"] - cenv["lactmo"]

df = df.rename(columns={
    "MomID.Unique": "MomIdUnique",
    "DadID.Unique": "DadIdUnique",
    "Stimulant_No_PerWeek": "stimulant",
})

# Replace these variables with their first value within a person.
df = df.sort_values(by=["ID", "datecomb"])
df["Wealth_Z"] = df.groupby("ID")["Wealth_Z"].transform("first")
df["Village"] = df.groupby("ID")["Village"].transform("first")

# Missing stimulant assumed to be no stimulant use
df["stimulant"] = df["stimulant"].fillna(0)

varls = [
    "ID", bp_var, "age_x", "Age", "Female",
    "Male", "BMI_cen", "HT_cen", "MomIdUnique",
    "Village"]

if adj >= 1:
    varls += ["Wealth_Z"]
    if sex == "female":
        varls += ["PregMo_Use_cen", "lactmo_Use_cen"]
    else:
        varls += ["stimulant"]

if adj >= 2:
    varls += ["Mom_Ht_Ave", "Mom_SBP_Z_Res_USE"]

if adj >= 3:
    varls += ["Dad_Ht_Ave", "Dad_SBP_Z_Res_USE"]

if sex == "female":
	df = df.loc[df.Sex == "Female", :]
elif sex == "male":
	df = df.loc[df.Sex == "Male", :]

log = open(os.path.join(bp_dir, d1, d2, "%s_med.log" % impvar), "w")

mn, cov, proj, vb, vx = do_pca(impvar, ndim, sex, bp_dir)

fn = os.path.join(bp_dir, d1, d2, "%s_%s_pca.txt" % (impvar, sex))
with open(fn, "w") as io:
    io.write(np.array2string(proj))
    io.write("\n===\n")
    io.write(np.array2string(np.sqrt(np.diag(cov))))

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

    bs_vars = ["HT_cen", "BMI_cen"]

    # Fit models for the mediators as outcomes
    rslt, fml = [], []
    scale = [None, None]
    for (j, bs_var) in enumerate(bs_vars):

        # Set up the formula
        f = bs_var + " ~ age_x + I(age_x**2) + I(age_x**3)"
        f += " + " + " + ".join(vb)
        f += " + C(Village)"
        if adj >= 1:
            f += " + Wealth_Z"
            if (sex == "female") and (bs_var == "BMI_cen"):
                f += " + PregMo_Use_cen + I(PregMo_Use_cen**2) + lactmo_Use_cen"
            if sex == "male":
                f += " + stimulant"
        if adj >= 2:
            f += " + Mom_Ht_Ave"
        if adj >= 3:
            f += " + Dad_Ht_Ave"
        fml.append(f)

        imp = make_imp(f)
        r = imp.fit(results_cb=lambda x: x)
        fv = np.dot(r.model.exog, r.params)
        scale[j] = np.mean([x.scale for x in r.results])
        rslt.append(r)

        out.write(r.summary().as_text())
        out.write("scale=%f\n" % scale[j])
        out.write("\n\n")

    # Standard deviations of the childhood trajectory
    s = np.sqrt(np.diag(cov))

    # Imputers for each mediator
    imp = [make_imp(f).imp for f in fml]

    #yv, mnd, mnv, pc = [], [], [], []
    ym, yv, ydm, ydv, pc = [], [], [], [], []
    for k in range(n_imp):

        # Get imputed design matrices for each mediator
        resid = [None, None]
        dz = []
        for j in 0,1: # HT, BMI
            imp[j].update()
            dx = imp[j].data
            dx = dx.groupby("ID").head(1)
            dx.loc[:, "age_x"] = (ref_age - 10) / 10
            dz.append(dx)
            xmat = patsy.dmatrix(rslt[j].model.data.design_info, data=dx,
                                 return_type='dataframe')
            resid[j] = dx[bs_vars[j]] - np.dot(xmat, rslt[j].params)

        # Get correlated errors for imputed HT and BMI
        r = np.corrcoef(resid[0], resid[1])[0, 1]
        n = dz[0].shape[0]
        ee = []
        for j in 0, 1: # HT, BMI
            e1 = np.random.normal(size=n)
            e2 = r*e1 + np.sqrt(1-r**2)*np.random.normal(size=n)
            ee.append([e1, e2])

        # Simulate data for both predictor variables under
        # two growth patterns.
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

        # Means and variances over subjects, within this imputation,
        # for HT (pat1, pat2) and BMI (pat1, pat2).
        ym.append([y.mean() for y in yy])
        yv.append([y.var() + scale[j//2] for (j,y) in enumerate(yy)])

        # pat2 - pat2 contrasts, averaged over subjects, within this
        # imputation.
        ydm.append([(yy[1] - yy[0]).mean(), (yy[3] - yy[2]).mean()])
        ydv.append([(yy[1] - yy[0]).var(), (yy[3] - yy[2]).var()])

    ym = np.asarray(ym)
    yv = np.asarray(yv)
    ydm = np.asarray(ydm)
    ydv = np.asarray(ydv)

    return ym, yv, ydm, ydv, (pc[0][0], pc[0][1]), patterns

ym, yv, ydm, ydv, pcl, patterns = do()
ymm = ym.mean(0)
ymv = yv.mean(0)
yvm = ym.var(0)
ydmm = ydm.mean(0)
ydvm = ydv.mean(0)
ydmv = ydm.var(0)

# Get the model coefficients
d = "mixed"
fn = os.path.join(bp_dir, d, "dim_%d" % ndim, "%s_adj%d_%s.txt" % (impvar, adj, sex))
coefs = {}
vcov = []
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
        line = line.replace("Village, Treatment", "Village,Treatment")
        line = line.replace(" ** ", "**")
        v = line.split()
        coefs[v[0]] = (float(v[1]), float(v[2]))

    # Get the parameter variance/covariance matrix
    for line in inf:
        vcov.append(line.split(","))
    vcov = vcov[-8:-2]
    for i in range(1, len(vcov)):
        for j in range(1, len(vcov[i])):
            vcov[i][j] = float(vcov[i][j])
    vcov = pd.DataFrame(vcov)
    vcov = vcov.set_index(0, drop=True)
    vcov.columns = vcov.iloc[0, :]
    vcov = vcov.iloc[1:, :]

out.write("Childhood %s %s versus %s\n" % (impvar, patterns[pat1][0], patterns[pat2][0]))
se = np.sqrt(ydmv[0] + ydvm[0])
out.write("Mean change in height at age %d:%10.3f cm (%.3f)\n" % (ref_age, ydmm[0], se))
se = np.sqrt(ydmv[1] + ydvm[1])
out.write("Mean change in BMI at age %d:   %10.3f kg/m^2 (%.3f)\n" % (ref_age, ydmm[1], se))

cf, cse = [], []
for s in ["HT_cen", "BMI_cen", "%s_pc0" % impvar, "HT_cen:%s_pc0" % impvar, "BMI_cen:%s_pc0" % impvar]:
    cf.append(coefs[s][0])
    cse.append(coefs[s][1])

def bs(v):
    va = np.dot(v, np.dot(vcov, v))
    return np.dot(cf, v), np.sqrt(va)

# Indirect effects through height
# [ymm[0], ymm[2], pcl[0], pcl[0] * ymm[0], pcl[0] * ymm[2]]
# [ymm[1], ymm[2], pcl[0], pcl[1] * ymm[1], pcl[0] * ymm[2]]
v = np.r_[ydmm[0], 0, 0, pcl[1] * ymm[1] - pcl[0] * ymm[0], 0]
inh, inh_se = bs(v)

# Indirect effects through BMI
# [ymm[0], ymm[2], pcl[0], pcl[0] * ymm[0], pcl[0] * ymm[2]]
# [ymm[0], ymm[3], pcl[0], pcl[0] * ymm[0], pcl[1] * ymm[3]]
v = np.r_[0, ydmm[1], 0, 0, pcl[1] * ymm[3] - pcl[0] * ymm[2]]
inb, inb_se = bs(v)

# Overall indirect effects
# [ymm[0], ymm[2], pcl[0], pcl[0] * ymm[0], pcl[0] * ymm[2]]
# [ymm[1], ymm[3], pcl[0], pcl[1] * ymm[1], pcl[1] * ymm[3]]
v = np.r_[ydmm[0], ydmm[1], 0, pcl[1] * ymm[1] - pcl[0] * ymm[0],
          pcl[1] * ymm[3] - pcl[0] * ymm[2]]
ino, ino_se = bs(v)

# Direct effects
# [ymm[0], ymm[2], pcl[0], pcl[0] * ymm[0], pcl[0] * ymm[2]]
# [ymm[0], ymm[2], pcl[1], pcl[0] * ymm[0], pcl[0] * ymm[2]]
v = np.r_[0, 0, pcl[1] - pcl[0], 0, 0]
de, de_se = bs(v)

# Total effects
# [ymm[0], ymm[2], pcl[0], pcl[0] * ymm[0], pcl[0] * ymm[2]]
# [ymm[1], ymm[3], pcl[1], pcl[1] * ymm[1], pcl[1] * ymm[3]]
v = np.r_[ydmm[0], ydmm[1], pcl[1] - pcl[0], pcl[1] * ymm[1] - pcl[0] * ymm[0],
          pcl[1] * ymm[3] - pcl[0] * ymm[2]]
te, te_se = bs(v)


a = inh
b = inh - 2*inh_se
c = inh + 2*inh_se
out.write("Indirect effect through height: %10.3f mm/Hg (%.3f,%.3f)\n" % (a, b, c))

a = inb
b = inb - 2*inb_se
c = inb + 2*inb_se
out.write("Indirect effect through BMI:    %10.3f mm/Hg (%.3f,%.3f)\n" % (a, b, c))

a = ino
b = ino - 2*ino_se
c = ino + 2*ino_se
out.write("Overall indirect effect:        %10.3f mm/Hg (%.3f,%.3f)\n" % (a, b, c))

a = de
b = de - 2*de_se
c = de + 2*de_se
out.write("Direct effect:                  %10.3f mm/Hg (%.3f,%.3f)\n" % (a, b, c))

a = te
b = te - 2*te_se
c = te + 2*te_se
out.write("Total effect:                   %10.3f mm/Hg (%.3f,%.3f)\n" % (a, b, c))

out.close()
