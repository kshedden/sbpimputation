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
if impvar not in ("HT", "WAZ", "HAZ", "BAZ", "WT"):
    msg = "Unknown imputation variable"
    sys.exit(msg)

df = pd.read_csv("/nfs/kshedden/Beverly_Strassmann/Cohort_2018.csv.gz")

# Selections
df = df.loc[df.Year > 2009]
df = df.loc[df.Bamako.isin([0, 1]), :]

df["Sex"] = df["Sex"].replace({0: "Female", 1: "Male"})
df["Female"] = df["Sex"].replace({"Female": 1, "Male": 0})
df["BMI_cen"] = df.BMI_15 - df.BMI_15.mean()
df["HT_cen"] = df.Ht_Ave - df.Ht_Ave.mean()
df["temp_cen"] = df.Temp - df.Temp.mean()
df["age_cen"] = df.Age_Yrs - df.Age_Yrs.mean()
df["age_x"] = (df.Age_Yrs - 10) / 10

df = df.rename(columns={
    "MomID.Unique": "MomIdUnique",
    "DadID.Unique": "DadIdUnique"
})

vars = [
    "ID", "SBMean23", "Bamako", "age_cen", "age_x", "Age_Yrs", "Female",
    "BMI_cen", "HT_cen", "Year", "lognummeas", "temp_cen", "School",
    "Wealth_Z", "MomIdUnique"
]

vx = [(impvar + "%d") % a for a in range(1, 11)]

low_thresh = {"HAZ": -2, "BAZ": -2, "WAZ": -2, "HT": 75, "WT": 8}[impvar]


# A class for producing imputed data sets
class mimi(object):

    ix = 0

    def update(self):
        di = pd.read_csv(
            os.path.join("imputed_data", "%s_imp_%d.csv" % (impvar, self.ix)))
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

        di["years_low"] = (dd < low_thresh).sum(1)
        di["min"] = dd.min(1)

        vb = ["ID"]
        vb += [(impvar + "_%s") % x for x in "ilq"]
        vb += [(impvar + "_d%s") % x for x in "ilq"]
        vb += ["years_low", "min"]
        di = di.loc[:, vb]

        dx = pd.merge(
            df.loc[:, vars], di, left_on="ID", right_on="ID", how="left")
        self.data = dx.dropna()


fml = "SBMean23 ~ Female*(age_x + I(age_x**2) + I(age_x**3)) + HT_cen + Female*BMI_cen + C(Year) + lognummeas + temp_cen + School + Wealth_Z + Bamako + "

fml_lin = fml + "%s_i + %s_l + %s_q" % (impvar, impvar, impvar)
fml_dlin = fml + "%s_di + %s_dl + %s_dq" % (impvar, impvar, impvar)
fml_min = fml + "min"
fml_yrlow = fml + "years_low"

out = open("imp_%s.txt" % impvar, "w")

# DEBUG
#vcf = {"Id": "0 + C(ID)", "Id*age": "0 + C(ID)*age_cen"}
vcf = {"Id": "0 + C(ID)"}


def model_kwds_fn(x):
    return {"groups": "MomIdUnique", "re_formula": "1", "vc_formula": vcf}


for fml in fml_dlin, fml_min, fml_yrlow, fml_lin:

    imp = MI(
        mimi(),
        sm.MixedLM,
        None,
        formula=fml,
        model_kwds_fn=model_kwds_fn,
        burn=0,
        nrep=20,
        skip=0)
    rslt = imp.fit()

    out.write("%d distinct subjects\n" % df.ID.unique().size)
    out.write(rslt.summary().as_text())

    # Save the coefficient trajectories, if present
    if fml in (fml_lin, fml_dlin):
        mp = rslt.params
        cm = rslt.cov_params()
        xn = rslt.model.exog_names
        p = len(xn)
        cm = cm[0:p, 0:p]
        xn = xn[0:p]
        x = "" if fml == fml_lin else "d"
        i = xn.index("%s_%si" % (impvar, x))
        l = xn.index("%s_%sl" % (impvar, x))
        q = xn.index("%s_%sq" % (impvar, x))
        ag = np.arange(1, 11)
        ax = np.linspace(-1, 1, 10)
        mx = np.zeros((len(ag), 3))
        for k, x in enumerate(ax):
            d = np.zeros(p)
            d[i] = 1
            d[l] = ax[k]
            d[q] = ax[k]**2
            mx[k, 0] = ag[k]
            mx[k, 1] = mp[i] + mp[l] * ax[k] + mp[q] * ax[k]**2
            mx[k, 2] = np.sqrt(np.dot(d, np.dot(cm, d)))
        mx = pd.DataFrame(mx, columns=["age", "coeff", "se"], index=ag)
        out.write(mx.to_string())

out.close()
