"""

Call the script with the name of the variable to be imputed, e.g.

> python impute_puberty.py log2T_use_Z
> python impute_puberty.py Breast_Stage_Use_Z
"""

import sys
sys.path.insert(
    0, "/afs/umich.edu/user/k/s/kshedden/statsmodels_fork/statsmodels")

import matplotlib
matplotlib.use('agg')
import numpy as np
import os
import pandas as pd
from data_tools import get_data, age_mean, age_sd
from statsmodels.regression.process_regression import ProcessMLE
import matplotlib.pyplot as plt
from config import *
from matplotlib.backends.backend_pdf import PdfPages
import statsmodels.api as sm

# Make sure that we are imputing a valid variable.
impvar = sys.argv[1]
if impvar not in ("log2T_use_Z", "Breast_Stage_Use_Z"):
    msg = "Unknown imputation variable"
    sys.exit(msg)

pdf = PdfPages("covmat_%s.pdf" % impvar)
outf = open("impute_%s_info.txt" % impvar, "w")

# Storage for results
dx = [None, None]
preg = [None, None]
rslt = [None, None]

# Carry values forward or backward within a time window
def xform(dz):
    ws = 90  # Days allowed to carry forward or backward
    imiss = pd.isnull(dz[impvar])
    if not imiss.any():
        return dz
    iobs = pd.notnull(dz[impvar])
    if not iobs.any():
        return dz
    dobs = dz.loc[iobs, :]
    for i in imiss.index[imiss]:
        td = (dz.loc[i, "datecomb"] - dobs.loc[iobs, "datecomb"]).dt.days
        ii = np.abs(td).idxmin()
        if np.abs(td[ii]) < ws:
            dz.loc[i, impvar] = dz.loc[ii, impvar]
    return dz

minage = 8
maxage = 20

# Additional variables to use when imputing puberty
# TODO: include Bamako, subscap for boys, men_age for girls
others = {0: ["HT", "WT", "BMI"],
          1: ["HT", "WT"]}

# Fit a Gaussian process model separately to females and males.
for female in 0, 1:

    if female == 1 and impvar in ("log2T_use_Z",):
        continue
    if female == 0 and impvar in ("Breast_Stage_Use_Z", "Menarche"):
        continue

    dx[female] = get_data(female, impvar, others=["datecomb"] + others[female])
    outf.write("Loaded %d x %d values\n" % tuple(dx[female].shape))
    outf.write("%d distinct people in initial data\n\n" % dx[female].ID.unique().size)

    dx[female] = dx[female].loc[dx[female].Age >= minage, :]
    outf.write("Retained %d x %d values at or above age %d\n" % (tuple(dx[female].shape) + (minage,)))
    outf.write("%d distinct people after requiring age at or above %d\n\n" % (dx[female].ID.unique().size, minage))

    # Not converted for some reason, stored as seconds, convert to days
    dx[female].datecomb = pd.to_datetime(dx[female].datecomb)

    dx[female] = dx[female].groupby("ID").apply(xform)

    # Drop if any variable other than impvar is missing
    ii = pd.notnull(dx[female].loc[:, others[female]]).all(1)
    dx[female] = dx[female].loc[ii, :]
    outf.write("Retained %d x %d values after dropping missing values\n" % tuple(dx[female].shape))
    outf.write("%d distinct people will have observed or imputed values\n\n" % dx[female].ID.unique().size)

    # Z-score age again
    age_mean = dx[female].Age.mean()
    age_sd = dx[female].Age.std()
    dx[female]["AgeZ"] = (dx[female].Age - age_mean) / age_sd

    ddm = dx[female].dropna()
    outf.write("%d x %d values used to fit model\n" % ddm.shape)
    outf.write("%d distinct people used to fit model\n\n" % ddm.ID.unique().size)

    # Fit the model -- note that degrees of freedom are fixed here.
    # TODO: are these good choices of the df values?
    # TODO: maybe this should be refit each imputation cycle with bootstrapped data
    preg[female] = ProcessMLE.from_formula(
        impvar + " ~ " + " + ".join(others[female]),
        scale_formula="AgeZ + I(AgeZ**2)",
        smooth_formula="1",
        noise_formula="1",
        time="Age",
        groups="ID",
        data=ddm)

    rslt[female] = preg[female].fit(verbose=True, maxiter=[2, 200])

    ages = np.linspace(minage, maxage, 100)
    dv = pd.DataFrame({"Age": ages, "AgeZ": (ages - age_mean) / age_sd})

    # Plot the mean function
    #minage = min(ages)
    #maxage = max(ages)
    #plt.clf()
    #plt.axes([0.1, 0.1, 0.7, 0.8])
    #for wtz in -1, 0, 1:
    #    dv = pd.DataFrame({"Age": ages, "AgeZ": (ages - age_mean) / age_sd, "WTZ": wtz})
    #    fv = rslt[female].predict(exog=dv)
    #    plt.plot(ages, fv, '-', lw=4, label="WTZ=%.0f" % wtz)
    #plt.grid(True)
    #plt.xlabel("Age", size=15)
    #plt.ylabel(impvar, size=15)
    #ha, lb = plt.gca().get_legend_handles_labels()
    #leg = plt.figlegend(ha, lb, "center right", handletextpad=0.0001)
    #leg.draw_frame(False)
    #pdf.savefig()

    # Plot the fitted covariance matrix
    mnpar = rslt[female].mean_params
    scpar = rslt[female].scale_params
    smpar = rslt[female].smooth_params
    cm = preg[female].covariance(ages, scpar, smpar, dv, dv)
    for k in 0, 1:
        # First plot covariance, then correlation
        plt.clf()
        plt.imshow(cm, interpolation='nearest', extent=[minage, maxage, maxage, minage])
        plt.xlabel("Age", size=15)
        plt.ylabel("Age", size=15)
        plt.title("%s (%s) %s" % (impvar, ["males", "females"][female], ["covariance", "correlation"][k]))
        plt.colorbar()
        pdf.savefig()
        s = np.sqrt(np.diag(cm))
        cm /= np.outer(s, s)

# Open files for storing all the imputed data.
out = []
for j in range(n_imp):
    out.append(
        open(os.path.join("imputed_data_puberty", "%s_imp_%d.csv" % (impvar, j)), "w"))

# Generate imputed data for females and males separately.
for female in 0, 1:

    if female == 1 and impvar in ("log2T_use_Z",):
        continue
    if female == 0 and impvar in ("Breast_Stage_Use_Z",):
        continue

    # Model parameters
    mnpar = rslt[female].mean_params
    scpar = rslt[female].scale_params
    smpar = rslt[female].smooth_params
    nopar = rslt[female].no_params

    # Remove the mean
    dx[female]["mean_impvar"] = rslt[female].predict(dx[female])
    dx[female]["resid"] = dx[female][impvar] - dx[female].mean_impvar

    # Loop over subjects
    first = True
    for idx, vx in dx[female].groupby("ID"):

        ix_obs = pd.notnull(vx[impvar]).values
        ix_miss = pd.isnull(vx[impvar]).values

        if ix_obs.all():
            # Just print the headers if needed
            for j in range(n_imp):
                out[j].write(vx.loc[ix_miss, ["ID", "Age", impvar]].to_csv(index=None, header=first))
            first = False
            continue

        print(idx)

        # Get the covariance matrix for all relevant ages.
        dv = pd.DataFrame({"Age": vx.Age, "AgeZ": vx.AgeZ})
        cm = preg[female].covariance(dv.Age.values, scpar, smpar, dv, dv)
        cm += np.exp(nopar[0])**2 * np.eye(cm.shape[0]) / 10
        cm00 = cm[ix_miss, :][:, ix_miss]
        cm01 = cm[ix_miss, :][:, ix_obs]
        cm11 = cm[ix_obs, :][:, ix_obs]

        # The mean for imputed values
        z = vx.mean_impvar[ix_miss] + np.dot(cm01, np.linalg.solve(cm11, vx.resid[ix_obs]))

        # The covariance matrix for imputed values
        va = cm00 - np.dot(cm01, np.linalg.solve(cm11, cm01.T))

        # The square root of the covariance matrix, warn if degenerate
        e, vm = np.linalg.eigh(va)
        e = e * (e >= 0)
        if np.any(e == 0):
            print("Covariance matrix is singular")
        vr = vm * np.sqrt(e)
        if np.iscomplex(vr[0,0]): 1/0

        # Generate the imputed values
        for j in range(n_imp):
            ze = z + np.dot(vr, np.random.normal(size=vr.shape[0]))
            vx.loc[ix_miss, impvar] = ze
            out[j].write(vx.loc[ix_miss, ["ID", "Age", impvar]].to_csv(index=None, header=first))
        first = False

# Close all the files
for j in range(n_imp):
    out[j].close()

pdf.close()
outf.close()