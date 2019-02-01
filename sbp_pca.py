"""
Fit models for outcome variables relative to exposures.

Pass the name of the outcome variable

> python sbp.py impvar [growth|nogrowth]
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
import json
from config import *

impvar = sys.argv[1]
if impvar not in allowed_controls:
    msg = "Unknown imputation variable"
    sys.exit(msg)

bp_var = "SBP_MEAN" if sys.argv[2] == "sbp" else "DBP_MEAN"
bp_dir = "sbp" if sys.argv[2] == "sbp" else "dbp"

mixed = sys.argv[3] == "mixed"

ndim = int(sys.argv[4])

growth = sys.argv[5] == "growth"

if growth and "Z" in impvar:
    sys.exit(0)

impvar_full = impvar + ("_growth" if growth else "")

# Save means before centering
cenv = {}

df["Male"] = df["Sex"].replace({"Female": 0, "Male": 1})
cenv["BMI"] = df.BMI_18.mean()
df["BMI_cen"] = df.BMI_18 - cenv["BMI"]
cenv["HT"] = df.Ht_Ave_18.mean()
df["HT_cen"] = df.Ht_Ave_18 - cenv["HT"]

df.Ht_Ave_18.mean()
cenv["temp"] = df.Temp.mean()
df["temp_cen"] = df.Temp - cenv["temp"]
cenv["age"] = df.Age.mean()
df["age_cen"] = df.Age - cenv["age"]
df["age_x"] = (df.Age - 10) / 10
df["School"] = df.School_Imputed

df["ID"] = df["ID"].astype(np.int)

df = df.rename(columns={
    "MomID.Unique": "MomIdUnique",
    "DadID.Unique": "DadIdUnique",
})

vars = [
    "ID", bp_var, "Bamako", "age_cen", "age_x", "Age", "Female",
    "Male", "BMI_cen", "HT_cen", "Year", "lognummeas", "temp_cen", "School",
    "Wealth_Z_2", "MomIdUnique", "log2T_use_Z", "Breast_Stage_Use_Z",
    "PregMo_Use_cen", "LactMo_Use_cen",
]

# These are missing for males
ii = (df.Female == 0)
df.loc[ii, "PregMo_Use"] = df.PregMo_Use[ii].fillna(0)
df.loc[ii, "LactMo_Use"] = df.LactMo_Use[ii].fillna(0)

cenv["LactMo"] = df.LactMo_Use.mean()
df["LactMo_Use_cen"] = df.LactMo_Use - cenv["LactMo"]
cenv["PregMo"] = df.PregMo_Use.mean()
df["PregMo_Use_cen"] = df.PregMo_Use - cenv["PregMo"]

# Save the means to a file
d1 = "mixed" if mixed else "gee"
d2 = "dim_%d" % ndim
fid = open(os.path.join(bp_dir, d1, d2, "%s_cen.json" % impvar_full), "w")
fid.write(json.dumps(cenv))
fid.close()

log = open(os.path.join(bp_dir, d1, d2, "%s.log" % impvar_full), "w")

# The variable names for the trajectory of childhood exposures
vx = [(impvar + "%d") % a for a in range(1, 11)]

# Second derivatives
p = len(vx)-2 if growth else len(vx)
if growth:
    p -= 1
fx = np.zeros((p-2, p))
for i in range(p-2):
    fx[i, i:i+3] = [1, -2, 1]

def get_growth(dd):
    dd = np.asarray(dd)
    dd = [np.correlate(y, np.r_[-0.3, -0.1, 0.1, 0.3], mode='valid') for y in dd]
    dd = np.asarray(dd)
    return dd

# Get the moments of the imputed data
mn, cov = 0, 0
for k in range(20):
    di = pd.read_csv(os.path.join("imputed_data", "%s_imp_%d.csv" % (impvar, k)))
    di["ID"] = di["ID"].astype(np.int)
    dd = di[vx]
    if growth:
        dd = get_growth(dd)
    mn += dd.mean(0)
    cov += np.cov(dd.T)
mn /= 20
cov /= 20

# Functional PCA
xsd = np.sqrt(np.diag(cov))
cor = cov / np.outer(xsd, xsd)
cm = cor - 0.2*np.dot(fx.T, fx)
b, proj = np.linalg.eig(cm)
ii = np.argsort(b)[::-1]
b = b[ii]
proj = proj[:, ii]
b = b[0:ndim]
proj = proj[:, 0:ndim]
proj = proj / xsd[:, None]
if np.any(b < 0):
    1/0

# Plot basis functions
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
d = "mixed" if mixed else "gee"
pdf = PdfPages(os.path.join(bp_dir, d, "dim_%d" % ndim, "%s_basis.pdf" % impvar_full))
plt.clf()
plt.title("%s basis functions" % impvar_full)
plt.grid(True)
for k in range(ndim):
    plt.plot(np.arange(1, proj.shape[0]+1), proj[:, k], '-', color='orange', lw=4)
plt.xlabel("Age", size=15)
pdf.savefig()
pdf.close()

vb = ["%s_pc%d" % (impvar, j) for j in range(proj.shape[1])]

# A class for producing imputed data sets
class mimi(object):

    ix = 0

    def update(self):

        # Create an array to hold the gridded values of impvar
        fn = os.path.join("imputed_data", "%s_imp_%d.csv" % (impvar, self.ix))
        di = pd.read_csv(fn)
        di["ID"] = di["ID"].astype(np.int)
        dd = di[vx]

        if growth:
            dd = get_growth(dd)

        dbs = np.dot(dd - mn, proj)
        dbs = np.concatenate((di.ID.values[:,None], dbs), axis=1)
        dbs = pd.DataFrame(dbs, columns=["ID"] + vb)
        dbs.ID = dbs.ID.astype(np.int)

        # Merge all variables created above with other cohort file variables
        dx = pd.merge(
            df.loc[:, vars], dbs, left_on="ID", right_on="ID", how="left")

        if self.ix == 0:
            log.write("%d distinct subjects in imputed data\n" % dbs.ID.unique().size)
            log.write("%d distinct subjects in cohort file data\n" % df.ID.unique().size)
            log.write("%d distinct subjects in merged file data\n" % dx.ID.unique().size)

        # Merge in imputed puberty variables
        for vn in "Breast_Stage_Use_Z", "log2T_use_Z":
            dd = pd.read_csv(os.path.join("imputed_data_puberty", "%s_imp_%d.csv" % (vn, self.ix)))
            dd = pd.merge(dx, dd, left_on=("ID", "Age"), right_on=("ID", "Age"), how="left")
            ix = pd.notnull(dd[vn + "_x"])
            iy = pd.notnull(dd[vn + "_y"])
            dd[vn] = np.nan
            dd.loc[ix, vn] = dd.loc[ix, vn + "_x"]
            dd.loc[iy, vn] = dd.loc[iy, vn + "_y"]
            dx = dd.drop([vn + "_x", vn + "_y"], axis=1)
            if self.ix == 0:
                log.write("%d distinct subjects after merging %s\n" %
                          (dx.ID.unique().size, vn))

        # Code testosterone for females, breast stage for males as zero
        dx.loc[dx.Female == 1, "log2T_use_Z"] = 0
        dx.loc[dx.Female == 0, "Breast_Stage_Use_Z"] = 0

        dx = dx.loc[pd.notnull(dx[bp_var]), :]
        self.data = dx.dropna()

        if self.ix == 0:
            log.write("%d distinct subjects in merged file data\n" % self.data.dropna().ID.unique().size)

        if self.ix == 0:
            dt = self.data[["ID", bp_var]].copy()
            nf = dt.groupby("ID")[bp_var].size()
            nf = pd.DataFrame(nf)
            nf = nf.reset_index()
            nf.columns = ["ID", "n_sbp"]
            method = {True: "mixed", False: "gee"}[mixed]
            nf.to_csv(os.path.join(bp_dir, method, "stats", "%s_n.csv" % impvar_full), index=False)

        self.ix += 1

fml = bp_var + " ~ Female*(age_x + I(age_x**2) + I(age_x**3) + lognummeas) + C(Year) + temp_cen + School + Wealth_Z_2 + Bamako + Female:Breast_Stage_Use_Z + Male:log2T_use_Z + PregMo_Use_cen + I(PregMo_Use_cen**2) + LactMo_Use_cen + I(LactMo_Use_cen**2)"

# Controls for body dimensions at time of SBP
fmx = fml + " + HT_cen + Female*BMI_cen + "

fmx += " + " + " + ".join(vb)
fmx += " + " + " + ".join(["HT_cen*" + x for x in vb])
fmx += " + " + " + ".join(["BMI_cen*" + x for x in vb])

d = "mixed" if mixed else "gee"
out = open(os.path.join(bp_dir, d, "dim_%d" % ndim, "%s.txt" % impvar_full), "w")

vcf = {"Id": "0 + C(ID)", "Id*age": "0 + C(ID)*age_cen"}

#vcf = {"Id": "0 + C(ID)"}

if mixed:
    def model_kwds_fn(x):
        return {"groups": "MomIdUnique", "data": x, "re_formula": "1", "vc_formula": vcf}
else:
    def model_kwds_fn(x):
        dep_data = x.ID.values[:, None]
        return {"groups": "MomIdUnique", "data": x, "dep_data": dep_data,
                "cov_struct": sm.cov_struct.Nested()}

if mixed:
    def fit_kwds_fn(x):
        return {"method": "lbfgs", "reml": False}
else:
    def fit_kwds_fn(x):
        return {"cov_type": "naive", "maxiter": 2}


imp = MI(
    mimi(),
    sm.MixedLM if mixed else sm.GEE,
    model_args_fn=None, #model_args_fn,
    formula=fmx,
    model_kwds_fn=model_kwds_fn,
    fit_kwds=fit_kwds_fn,
    burn=0,
    nrep=20,
    skip=0)

rslt = imp.fit(results_cb=lambda x: x)

if mixed:
    ic = [x.llf for x in rslt.results]
else:
    ic = [x.qic(scale=93)[1] for x in rslt.results]

sca = [x.scale for x in rslt.results]

out.write("%s\n" % impvar_full)
out.write("%d distinct subjects\n" % df.ID.unique().size)
out.write("mean IC %f\n" % np.mean(ic))
out.write("mean scale %f\n" % np.mean(sca))
out.write(rslt.summary().as_text())

# Save the coefficient trajectories
vb1 = ["HT_cen", "BMI_cen", "Female:BMI_cen"] + vb + ["HT_cen:%s" % y for y in vb] + ["BMI_cen:%s" % y for y in vb]
mp = rslt.params
cm = rslt.cov_params()
xn = rslt.model.exog_names
mp = pd.Series(mp[0:len(xn)], xn)
mp = mp[vb1]
cm = pd.DataFrame(cm[0:len(xn), 0:len(xn)], index=xn, columns=xn)
cm = cm.loc[vb1, vb1]

# For growth, the age doesn't start at 1, might be better to change that
# here.
out.write("BEGIN-PROJECTION\n")
for k in range(proj.shape[0]):
    out.write(",".join([str(x) for x in proj[k, :]]) + "\n")
out.write("--\n")
for v in np.diag(cov):
    out.write("%f\n" % v)
out.write("--\n")
out.write(mp.to_csv())
out.write("--\n")
out.write(cm.to_csv())
out.write("END-PROJECTION\n\n")

if not mixed:
    out.write("%s\n" % rslt.results[0].cov_struct.summary())

out.close()
log.close()