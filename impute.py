"""
Generate imputed data for one stature variable.  Impute to a grid of ages, 1, 2, ... maxage,
where maxage is defined in config.py.

Call the script with the name of the variable to be imputed, e.g.

> python impute.py HT
> python impute.py BAZ
"""

import sys
sys.path.insert(
    0, "/afs/umich.edu/user/k/s/kshedden/statsmodels_fork/statsmodels")

import matplotlib
matplotlib.use('agg')
import numpy as np
import os
import pandas as pd
from data_tools import get_data
from statsmodels.regression.process_reg import ProcessRegression
import matplotlib.pyplot as plt
from config import *
from matplotlib.backends.backend_pdf import PdfPages

# Make sure that we are imputing a valid variable.
impvar = sys.argv[1]
if impvar not in allowed_controls:
    msg = "Unknown imputation variable"
    sys.exit(msg)

pdf = PdfPages("covmat_%s.pdf" % impvar)

# Ages to impute
imp_ages = np.arange(1, maxage + 1)

# Storage for results
dx = [None, None]
preg = [None, None]
rslt = [None, None]

# Fit a Gaussian process model separately to females and males.
for female in 0, 1:

    dx[female] = get_data(female, impvar)
    dx[female] = dx[female].loc[dx[female].Age <= maxage+1, :]
    dx[female] = dx[female].dropna()

    # Fit the model -- note that degrees of freedom are fixed here.
    # TODO: are these good choices of the df values?
    # TODO: maybe this should be refit each imputation cycle with bootstrapped data
    preg[female] = ProcessRegression.from_formula(
        "%s ~ bs(Age, 4)" % impvar,
        scale_formula="bs(Age, 4)",
        smooth_formula="bs(Age, 4)",
        time="Age",
        groups="ID",
        data=dx[female])

    rslt[female] = preg[female].fit()

    # Plot the fitted covariance matrix
    ages = np.linspace(1, maxage, 100)
    dv = pd.DataFrame({"Age": ages})
    mnpar = rslt[female].mean_params
    scpar = rslt[female].scale_params
    smpar = rslt[female].smooth_params
    cm = preg[female].covariance(ages, scpar, smpar, dv, dv)
    for k in 0, 1:
        # First plot covariance, then correlation
        plt.clf()
        plt.imshow(cm, interpolation='nearest', extent=[1, maxage, maxage, 1])
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
        open(os.path.join("imputed_data", "%s_imp_%d.csv" % (impvar, j)), "w"))
    out[j].write("ID," + ",".join(["%s%d" % (impvar, a) for a in imp_ages]))
    out[j].write(",")
    out[j].write(",".join(["Mean%d" % a for a in imp_ages]))
    out[j].write("\n")

# Generate imputed data for females and males separately.
for female in False, True:

    # Model parameters
    mnpar = rslt[female].mean_params
    scpar = rslt[female].scale_params
    smpar = rslt[female].smooth_params

    # Remove the mean
    exog = pd.DataFrame({"Age": imp_ages})
    mean_impvar = rslt[female].predict(exog)
    dx[female]["resid"] = dx[female][impvar] - rslt[female].predict()

    # Loop over subjects
    for idx, v in dx[female].groupby("ID"):

        # Create an age variable for all the ages we are concerned
        # with for this subject -- first the observed ages, then the
        # ages to be imputed.
        ages = v.Age.copy().values
        ages = np.concatenate((ages, imp_ages))
        dv = pd.DataFrame({"Age": ages})

        # Get the covariance matrix for all relevant ages.
        cm = preg[female].covariance(ages, scpar, smpar, dv, dv)
        p = dv.shape[0] - len(imp_ages)
        cm00 = cm[p:, p:]
        cm01 = cm[p:, 0:p]
        cm11 = cm[0:p, 0:p]

        # The mean for imputed values
        z = mean_impvar + np.dot(cm01, np.linalg.solve(cm11, v.resid))

        # The covariance matrix for imputed values
        va = cm00 - np.dot(cm01, np.linalg.solve(cm11, cm01.T))

        # The square root of the covariance matrix, warn if degenerate
        e, v = np.linalg.eig(va)
        e = e * (e >= 0)
        if np.any(e == 0):
            print("Covariance matrix is singular")
        vr = v * np.sqrt(e)

        # Generate the imputed values
        for j in range(n_imp):
            ze = z + np.dot(vr, np.random.normal(size=vr.shape[0]))
            out[j].write("%d," % idx)
            out[j].write(",".join(["%.3f" % u for u in ze]))
            out[j].write(",")
            out[j].write(",".join(["%.3f" % u for u in mean_impvar]))
            out[j].write("\n")

# Close all the files
for j in range(n_imp):
    out[j].close()

pdf.close()