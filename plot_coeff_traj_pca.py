"""
Generate plots of the time-varying coefficients for stature variables
as functions of age.
"""

import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import os
import matplotlib.pyplot as plt
from config import *
import sys
from io import StringIO

if len(sys.argv) != 3:
    raise ValueError("plot_coeff_traj: incorrect number of arguments")

bp_dir = "sbp" if sys.argv[1] == "sbp" else "dbp"

di = "mixed"

ndim = int(sys.argv[2])

pdf = PdfPages(os.path.join(bp_dir, di, "dim_%d" % ndim, "coeff_traj_pca.pdf"))

ages = np.arange(1, maxage + 1)

for vn in allowed_controls:
    for controlcbs in False, True:
        for dadbp in False, True:

            fn = os.path.join(bp_dir, di, "dim_%d" % ndim, "%s_nocontrolcbs_nodadbp.txt" % vn)
            if controlcbs:
                fn = fn.replace("nocontrolcbs", "controlcbs")
            if dadbp:
                fn = fn.replace("nodadbp", "dadbp")
            f = open(fn)

            cf = []
            for line in f:
                if line.startswith("BEGIN-PROJECTION"):
                    break

            cx = []
            for line in f:
                if line.startswith("--") or line.startswith("END-PROJECTION"):
                    cf.append(cx)
                    cx = []
                    continue
                if line.startswith("END-PROJECTION"):
                    break
                cx.append(line.rstrip())

            proj = pd.read_csv(StringIO("\n".join(cf[0])), header=None)
            vx = pd.read_csv(StringIO("\n".join(cf[1])), header=None)

            params = "\n".join(cf[2])
            params = pd.read_csv(StringIO(params), header=None, index_col=0)
            params = params[1]

            cov = "\n".join(cf[3])
            cov = pd.read_csv(StringIO(cov), index_col=0)

            age = np.arange(1, len(proj)+1, dtype=np.float64)

            if ndim == 1:
                vn1 = vn + "_pc0"
                traj = params[vn1] * proj
                c = cov.loc[vn1, vn1]
                va = np.diag(c * np.dot(proj, proj.T))
                se = np.sqrt(va)
                traj = traj.values.ravel()
            else:
                vn1 = vn + "_pc0"
                vn2 = vn + "_pc1"
                vn12 = [vn1, vn2]
                pa = params[vn12]
                traj = np.dot(proj, pa)
                c = cov.loc[vn12, vn12]
                va = np.diag(np.dot(proj, np.dot(c, proj.T)))
                se = np.sqrt(va)

            plt.clf()
            plt.axes([0.15, 0.11, 0.8, 0.8])
            plt.grid(True)
            plt.plot(age, traj, '-', color='black', lw=4)
            plt.fill_between(
                age,
                traj.flat - 2 * se,
                traj.flat + 2 * se,
                color='grey',
                alpha=0.5)
            plt.plot([1, maxage], [0, 0], '-', color='black', lw=3)
            plt.xlabel("Age in years", size=15)
            plt.ylabel("Coefficient for %s" % vn, size=15)

            title = "No control for height and BMI at time of SBP"
            if controlcbs:
                title = "Control for height and BMI at time of SBP"

            if dadbp:
                title += ", with dad BP"
            else:
                title += ", no dad BP"

            plt.title(title)
            pdf.savefig()

pdf.close()
