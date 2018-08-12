"""
Fit models for outcome variables relative to exposures.

Pass the name of the outcome variable

> python sbp.py impvar version

Versions define different model structures

1:
Control for height and BMI * gender

2:
Control for current impvar instead of height and BMI * gender
"""

import sys
sys.path.insert(
    0, "/afs/umich.edu/user/k/s/kshedden/statsmodels_fork/statsmodels")

import pandas as pd
import statsmodels.api as sm
from statsmodels.imputation.bayes_mi import MI
import os
import numpy as np
from config import *

impvar = sys.argv[1]
if impvar not in allowed_controls:
    msg = "Unknown imputation variable"
    sys.exit(msg)

version = int(sys.argv[2])

print("!!! %s %d\n" % (impvar, version))

df = pd.read_csv("/nfs/kshedden/Beverly_Strassmann/Cohort_2018_080118.csv.gz")

# Selections
df = df.loc[df.Bamako.isin([0, 1]), :]

df["Sex"] = df["Sex"].replace({0: "Female", 1: "Male"})
df["Female"] = df["Sex"].replace({"Female": 1, "Male": 0})
df["Male"] = df["Sex"].replace({"Female": 0, "Male": 1})
df["BMI_cen"] = df.BMI_18 - df.BMI_18.mean()
df["HT_cen"] = df.Ht_Ave_18 - df.Ht_Ave_18.mean()
df["temp_cen"] = df.Temp - df.Temp.mean()
df["Age"] = np.around(df["Age_Yrs"], 3)
df["age_cen"] = df.Age - df.Age.mean()
df["age_x"] = (df.Age - 10) / 10
df["HT"] = df.Ht_Ave_18
df["HAZ"] = df.HAZ_18
df["BAZ"] = df.BAZ_18
df["WAZ"] = df.WAZ_18
df["School"] = df.School_Imputed

df["log2T_use"] = np.log(df.T_use) / np.log(2)
df["log2T_use_Z"] = (df.log2T_use - df.log2T_use.mean()) / df.log2T_use.std()

df["Breast_Stage_Z"] = (df.Breast_Stage - df.Breast_Stage.mean()) / df.Breast_Stage.std()

df["ID"] = df["ID"].astype(np.int)

df = df.rename(columns={
    "MomID.Unique": "MomIdUnique",
    "DadID.Unique": "DadIdUnique",
})

vars = [
    "ID", "SBP_MEAN", "Bamako", "age_cen", "age_x", "Age", "Female",
    "Male", "BMI_cen", "HT_cen", "Year", "lognummeas", "temp_cen", "School",
    "Wealth_Z", "MomIdUnique", "log2T_use_Z", "Breast_Stage_Z",
]

if version == 2:
    vars.append(impvar)

vx = [(impvar + "%d") % a for a in range(1, 11)]

low_thresh = {"HAZ": -2, "BAZ": -2, "WAZ": -2, "HT": 75, "WT": 8}[impvar]


# A class for producing imputed data sets
class mimi(object):

    ix = 0

    def update(self):
        di = pd.read_csv(
            os.path.join("imputed_data", "%s_imp_%d.csv" % (impvar, self.ix)))
        di["ID"] = di["ID"].astype(np.int)
        dd = di[vx]
        q = dd.shape[1]
        di["%s_i" % impvar] = np.dot(dd, np.ones(q))
        di["%s_l" % impvar] = np.dot(dd, np.linspace(-1, 1, q))
        di["%s_q" % impvar] = np.dot(dd, np.linspace(-1, 1, q)**2)

        # Derivatives
        d1 = dd.diff(1, axis=1).iloc[:, 1:]
        q = d1.shape[1]
        di["%s_di" % impvar] = np.dot(d1, np.ones(q))
        di["%s_dl" % impvar] = np.dot(d1, np.linspace(-1, 1, q))
        di["%s_dq" % impvar] = np.dot(d1, np.linspace(-1, 1, q)**2)

        # Derivatives of logs (relative growth)
        if "Z" not in impvar:
            d1 = np.log(dd).diff(1, axis=1).iloc[:, 1:]
            q = d1.shape[1]
            di["%s_ri" % impvar] = np.dot(d1, np.ones(q))
            di["%s_rl" % impvar] = np.dot(d1, np.linspace(-1, 1, q))
            di["%s_rq" % impvar] = np.dot(d1, np.linspace(-1, 1, q)**2)

        di["years_low"] = (dd < low_thresh).sum(1)
        di["min"] = dd.min(1)

        vb = ["ID", "%s10" % impvar]
        vb += [(impvar + "_%s") % x for x in "ilq"]
        vb += [(impvar + "_d%s") % x for x in "ilq"]
        if "Z" not in impvar:
            vb += [(impvar + "_r%s") % x for x in "ilq"]
        vb += ["years_low", "min"]
        di = di.loc[:, vb]

        dx = pd.merge(
            df.loc[:, vars], di, left_on="ID", right_on="ID", how="left")

        for vn in "Breast_Stage_Z", "log2T_use_Z":
            dd = pd.read_csv(os.path.join("imputed_data_puberty", "%s_imp_%d.csv" % (vn, self.ix)))
            dd = pd.merge(dx, dd, left_on=("ID", "Age"), right_on=("ID", "Age"), how="left")
            dd[vn] = np.nan
            ix = pd.notnull(dd[vn + "_x"])
            iy = pd.notnull(dd[vn + "_y"])
            dd.loc[ix, vn] = dd.loc[ix, vn + "_x"]
            dd.loc[iy, vn] = dd.loc[iy, vn + "_y"]
            dx = dd.drop([vn + "_x", vn + "_y"], axis=1)

        dx.loc[dx.Female == 1, "log2T_use_Z"] = 0
        dx.loc[dx.Female == 0, "Breast_Stage_Z"] = 0

        self.data = dx.dropna()
        self.ix += 1

fml = "SBP_MEAN ~ Female*(age_x + I(age_x**2) + I(age_x**3)) + C(Year) + lognummeas + temp_cen + School + Wealth_Z + Bamako + Female:Breast_Stage_Z + Male:log2T_use_Z + "

if version == 1:
    fml += "HT_cen + Female*BMI_cen + "
elif version == 2:
    fml += "%s + " % impvar

fml_lin = fml + "%s_i + %s_l + %s_q" % (impvar, impvar, impvar)
fml_dlin = fml + "%s_di + %s_dl + %s_dq" % (impvar, impvar, impvar)
fml_rlin = fml + "%s_ri + %s_rl + %s_rq" % (impvar, impvar, impvar)
fml_min = fml + ("%s10 + min" % impvar)
fml_yrlow = fml + ("%s10 + years_low" % impvar)

out = open("imp_%s_%d.txt" % (impvar, version), "w")

vcf = {"Id": "0 + C(ID)", "Id*age": "0 + C(ID)*age_cen"}

#vcf = {"Id": "0 + C(ID)"}


def model_kwds_fn(x):
    return {"groups": "MomIdUnique", "re_formula": "1", "vc_formula": vcf}

def fit_kwds_fn(x):
    return {"method": "lbfgs"}

for fml in fml_dlin, fml_rlin, fml_min, fml_yrlow, fml_lin:

    if "Z" in impvar and fml == fml_rlin:
        continue

    imp = MI(
        mimi(),
        sm.MixedLM,
        None,
        formula=fml,
        model_kwds_fn=model_kwds_fn,
        fit_kwds=fit_kwds_fn,
        burn=0,
        nrep=20,
        skip=0)

    rslt = imp.fit()

    out.write("%s\n" % impvar)
    out.write("%d distinct subjects\n" % df.ID.unique().size)
    out.write(rslt.summary().as_text())

    # Save the coefficient trajectories, if present
    if fml in (fml_lin, fml_dlin, fml_rlin):
        mp = rslt.params
        cm = rslt.cov_params()
        xn = rslt.model.exog_names
        p = len(xn)
        cm = cm[0:p, 0:p]
        xn = xn[0:p]
        x = ""
        if fml == fml_dlin:
            x = "d"
        elif fml == fml_rlin:
            x = "r"
        i = xn.index("%s_%si" % (impvar, x))
        l = xn.index("%s_%sl" % (impvar, x))
        q = xn.index("%s_%sq" % (impvar, x))
        ax = np.linspace(-1, 1, 10)
        mx = np.zeros((10, 3))
        for k, x in enumerate(ax):
            d = np.zeros(p)
            d[i] = 1
            d[l] = ax[k]
            d[q] = ax[k]**2
            mx[k, 0] = k + 1
            mx[k, 1] = mp[i] + mp[l] * ax[k] + mp[q] * ax[k]**2
            mx[k, 2] = np.sqrt(np.dot(d, np.dot(cm, d)))
        mx = pd.DataFrame(mx, columns=["age", "coeff", "se"], index=range(1, 11))
        out.write("BEGIN-TRAJECTORY\n")
        out.write(mx.to_string(index=False))
        out.write("\nEND-TRAJECTORY\n")

out.close()
