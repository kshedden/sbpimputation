import sys
sys.path.insert(0, "/afs/umich.edu/user/k/s/kshedden/statsmodels_fork/statsmodels")

import pandas as pd
import numpy as np
from statsmodels.stats.correlation_tools import (kernel_covariance, corr_nearest)
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import statsmodels.api as sm

pdf = PdfPages("sbp_haz_kernel.pdf")

from sbp_data import get_data

for female in False, True:

    df = get_data(female)

    # Remove the mean structure from HAZ
    model = sm.OLS.from_formula("HAZ ~ bs(Age, 5)", data=df)
    result = model.fit()
    df["HAZ_resid"] = result.resid

    # Remove the mean structure from SBP
    model = sm.OLS.from_formula("SBMean23 ~ bs(Age, 5)", data=df)
    result = model.fit()
    df["SBMean23_resid"] = result.resid

    age = df.Age.values

    # Estimate the covariances and cross covariances among HAZ
    # and SBP
    cv = kernel_covariance(exog=df[["HAZ_resid", "SBMean23_resid"]],
                           loc=age, groups=df.ID, bw=2)

    # Evaluate the estimated covariances on a grid of ages
    agx = np.linspace(age.min(), age.max(), 40)
    cmhh = np.zeros((len(agx), len(agx)))
    cmbb = np.zeros((len(agx), len(agx)))
    cmhb = np.zeros((len(agx), len(agx)))
    print("Processing %s..." % {False: "males", True: "females"}[female])
    for i, a1 in enumerate(agx):
        print("%d/%d" % (i, len(agx)))
        for j, a2 in enumerate(agx):
            c = cv(a1, a2)
            cmhh[i, j] = c[0, 0]
            cmbb[i, j] = c[1, 1]
            cmhb[i, j] = c[0, 1]

    # Project to PSD
    if False:
        p = cmhh.shape[0]
        mat = np.zeros((2*p, 2*p))
        mat[0:p, 0:p] = cmhh
        mat[0:p, p:] = cmhb
        mat[p:, 0:p] = cmhb.T
        mat[p:, p:] = cmbb
        d = np.sqrt(np.diag(mat))
        mat *= np.outer(1/d, 1/d)
        mat = corr_clipped(mat) # corr_nearest
        mat *= np.outer(d, d)
        a,b = np.linalg.eig(mat)
        cmbb = mat[0:p, 0:p]
        cmhb = mat[0:p, p:]
        cmbb = mat[p:, p:]

    # Obtain correlation matrices
    d1 = np.sqrt(np.diag(cmhh))
    cmhh_c = cmhh * np.outer(1/d1, 1/d1)
    d2 = np.sqrt(np.diag(cmbb))
    cmbb_c = cmbb * np.outer(1/d2, 1/d2)
    cmhb_c = cmhb * np.outer(1/d1, 1/d2)

    for jx in 0,1:
        mat = [cmhh, cmhh_c][jx]
        ti = ["covariance", "correlation"][jx]
        plt.clf()
        plt.title(["Male ", "Female "][female] + ti)
        kw = {}
        if jx == 1:
            kw = {"vmin": 0.8, "vmax": 1}
        plt.imshow(mat, interpolation='nearest',
                   extent=[age.min(), age.max(), age.max(), age.min()],
                   **kw)
        plt.xlabel("HAZ age", size=16)
        plt.ylabel("HAZ age", size=16)
        plt.colorbar()
        pdf.savefig()

    for jx in 0,1:
        mat = [cmbb, cmbb_c][jx]
        ti = ["covariance", "correlation"][jx]
        plt.clf()
        plt.title(["Male ", "Female "][female] + ti)
        kw = {}
        if jx == 1:
            kw = {"vmin": 0.5, "vmax": 1}
        plt.imshow(mat, interpolation='nearest',
                   extent=[age.min(), age.max(), age.max(), age.min()],
                   **kw)
        plt.xlabel("SBP age", size=16)
        plt.ylabel("SBP age", size=16)
        plt.colorbar()
        pdf.savefig()

    for jx in 0,1:
        mat = [cmhb, cmhb_c][jx]
        ti = ["cross covariance", "cross correlation"][jx]
        plt.clf()
        plt.title(["Male ", "Female "][female] + ti)
        kw = {}
        if jx == 1:
            kw = {"vmin": 0, "vmax": 1}
        plt.imshow(mat, interpolation='nearest',
                   extent=[age.min(), age.max(), age.max(), age.min()],
                   **kw)
        plt.xlabel("SBP age", size=16)
        plt.ylabel("HAZ age", size=16)
        plt.colorbar()
        pdf.savefig()

pdf.close()