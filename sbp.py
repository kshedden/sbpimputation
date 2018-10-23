"""
Fit models for outcome variables relative to exposures.

Pass the name of the outcome variable

> python sbp.py impvar version

Versions define different model structures

1:
Control for height and BMI * gender

2:
Control for current impvar instead of height and BMI * gender

heat map !!!

T_date_use : saliva date

include reproductive status (reprostat)

go back to imputing testosterone

sports, smoking, minutes walked
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
from config import *

impvar = sys.argv[1]
if impvar not in allowed_controls:
    msg = "Unknown imputation variable"
    sys.exit(msg)

version = int(sys.argv[2])

print("!!! %s %d\n" % (impvar, version))

df["Male"] = df["Sex"].replace({"Female": 0, "Male": 1})
df["BMI_cen"] = df.BMI_18 - df.BMI_18.mean()
df["HT_cen"] = df.Ht_Ave_18 - df.Ht_Ave_18.mean()
df["temp_cen"] = df.Temp - df.Temp.mean()
df["age_cen"] = df.Age - df.Age.mean()
df["age_x"] = (df.Age - 10) / 10
df["School"] = df.School_Imputed


df["ID"] = df["ID"].astype(np.int)

df = df.rename(columns={
    "MomID.Unique": "MomIdUnique",
    "DadID.Unique": "DadIdUnique",
})

vars = [
    "ID", "SBP_MEAN", "Bamako", "age_cen", "age_x", "Age", "Female",
    "Male", "BMI_cen", "HT_cen", "Year", "lognummeas", "temp_cen", "School",
    "Wealth_Z", "MomIdUnique", "log2T_use_Z", "Breast_Stage_Use_Z",
    "Pregnant", "Lactating",
]

if version == 2:
    vars.append(impvar)

vx = [(impvar + "%d") % a for a in range(1, 11)]

low_thresh = {"HAZ": -2, "BAZ": -2, "WAZ": -2, "HT": 75, "WT": 8, "BMI": 15}[impvar]

import patsy
bs = patsy.dmatrix("0 + bs(a, 5)", pd.DataFrame({"a": np.linspace(-1, 1, 10)}), return_type='dataframe')
bs = np.asarray(bs)
bs[:, 1:] -= bs[:, 1:].mean(0)
bs /= np.sqrt((bs * bs).sum(0))


# Get the moments of the projected data
mnx, mn, cov = 0, 0, 0
for k in range(20):
    di = pd.read_csv(os.path.join("imputed_data", "%s_imp_%d.csv" % (impvar, k)))
    di["ID"] = di["ID"].astype(np.int)
    dd = di[vx]
    mnx += dd.mean(0)
    dbs = np.dot(dd, bs)
    mn += dbs.mean(0)
    cov += np.cov(dbs.T)
mn /= 20
mnx /= 20
cov /= 20
b, proj = np.linalg.eig(cov)


def xfeat(z):
    #z = [np.arctan(z+c) for c in (0,)]
    #z = np.concatenate(z, axis=1)
    return z

vbl = 5
vb = ["%s_x%d" % (impvar, j) for j in range(vbl)]

# A class for producing imputed data sets
class mimi(object):

    ix = 0

    def update(self):

        # Load the gridded values of impvar
        di = pd.read_csv(
            os.path.join("imputed_data", "%s_imp_%d.csv" % (impvar, self.ix)))
        di["ID"] = di["ID"].astype(np.int)
        dd = di[vx]

        dbs = np.dot(dd - mnx, bs)
        dbs = np.dot(dbs, proj)
        dbs = xfeat(dbs)
        dbs = np.concatenate((di.ID.values[:,None], dbs), axis=1)

        dbs = pd.DataFrame(dbs, columns=["ID"] + vb)
        dbs.ID =dbs.ID.astype(np.int)

        # Merge all variables created above with other cohort file variables
        dx = pd.merge(
            df.loc[:, vars], dbs, left_on="ID", right_on="ID", how="left")

        # Merge in imputed puberty variables
        for vn in "Breast_Stage_Use_Z", "log2T_use_Z":
            dd = pd.read_csv(os.path.join("imputed_data_puberty", "%s_imp_%d.csv" % (vn, self.ix)))
            dd = pd.merge(dx, dd, left_on=("ID", "Age"), right_on=("ID", "Age"), how="left")
            ix = pd.notnull(dd[vn + "_x"])
            iy = pd.notnull(dd[vn + "_y"])
            #assert(!(ix & iy).any())
            dd[vn] = np.nan
            dd.loc[ix, vn] = dd.loc[ix, vn + "_x"]
            dd.loc[iy, vn] = dd.loc[iy, vn + "_y"]
            dx = dd.drop([vn + "_x", vn + "_y"], axis=1)

        # Code testosterone for females, breast stage for males as zero
        dx.loc[dx.Female == 1, "log2T_use_Z"] = 0
        dx.loc[dx.Female == 0, "Breast_Stage_Use_Z"] = 0

        self.data = dx.dropna()
        self.ix += 1

fml = "SBP_MEAN ~ Female*(age_x + I(age_x**2) + I(age_x**3)) + C(Year) + lognummeas + temp_cen + School + Wealth_Z + Bamako + Female:Breast_Stage_Use_Z + Male:log2T_use_Z + Pregnant + Lactating + "

if version == 1:
    fml += "HT_cen + Female*BMI_cen + "
elif version == 2:
    fml += "%s + " % impvar

out = open("imp_%s_%d_x.txt" % (impvar, version), "w")

#vcf = {"Id": "0 + C(ID)", "Id*age": "0 + C(ID)*age_cen"}

vcf = {"Id": "0 + C(ID)"}


def model_kwds_fn(x):
    return {"groups": "MomIdUnique", "re_formula": "1", "vc_formula": vcf}

def fit_kwds_fn(x):
    return {"method": "lbfgs"}

fml += " + ".join(vb)

imp = MI(
    mimi(),
    sm.OLS, #DEBUG sm.MixedLM,
    None,
    formula=fml,
    #DEBUG model_kwds_fn=model_kwds_fn,
    #DEBUG fit_kwds=fit_kwds_fn,
    burn=0,
    nrep=20,
    skip=0)

rslt = imp.fit()

out.write("%s\n" % impvar)
out.write("%d distinct subjects\n" % df.ID.unique().size)
out.write(rslt.summary().as_text())

tsts = [np.linspace(-10, 0, len(mnx)), np.linspace(0, -10, len(mnx))]

pr = []
for tst in tsts:

    ts = np.dot(np.dot(tst, bs), proj)

    mi = mimi()
    pr0 = 0
    for k in range(20):

        mi.update()
        da = mi.data
        dx = da.iloc[0:10, :].copy()
        m = dx.mean(0)
        for i in range(10):
            dx.iloc[i, :] = m
        dx = dx.reset_index()
        dx.Year = 2012
        dx0 = dx.copy()
        dx1 = dx.copy()

        z = np.outer(np.linspace(0, 1, 10), ts)
        dx1.loc[:, vb] = xfeat(z)

        dx0.loc[:, vb] = xfeat(0* z)

        pr0 += rslt.predict(exog=dx1) - rslt.predict(exog=dx0)

    pr0 /= 20
    pr.append(pr0)
1/0

# Save the coefficient trajectories, if present
mp = rslt.params
cm = rslt.cov_params()
xn = rslt.model.exog_names
p = len(xn)
cm = cm[0:p, 0:p]
xn = xn[0:p]


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
out.write("\nEND-TRAJECTORY\n\n")

out.close()
