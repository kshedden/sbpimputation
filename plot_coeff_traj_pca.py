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

if len(sys.argv) != 4:
    raise ValueError("plot_coeff_traj: incorrect number of arguments")

bp_dir = "sbp" if sys.argv[1] == "sbp" else "dbp"

mixed = sys.argv[2] == "mixed"
di = "mixed" if mixed else "gee"

ndim = int(sys.argv[3])

pdf = PdfPages(os.path.join(bp_dir, di, "dim_%d" % ndim, "coeff_traj_pca.pdf"))

ages = np.arange(1, maxage + 1)

for vn in allowed_controls:
    for growth in False, True:

        vn_full = vn + ("_growth" if growth else "")

        if growth and "Z" in vn:
            continue

        fn = os.path.join(bp_dir, di, "dim_%d" % ndim, "%s.txt" % vn_full)
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

        cf[0] = [float(x) for x in cf[0]]
        proj = np.asarray(cf[0])

        cf[1] = [float(x) for x in cf[1]]
        vx = np.asarray(cf[1])

        from io import StringIO

        params = "\n".join(cf[2])
        params = pd.read_csv(StringIO(params), header=None, index_col=0)
        params = params[1]

        cov = "\n".join(cf[3])
        cov = pd.read_csv(StringIO(cov), index_col=0)

        age = np.arange(1, len(proj)+1, dtype=np.float64)
        if growth:
            age += 1.5

        vn1 = vn + "_pc0"
        traj = params[vn1] * proj
        se = np.sqrt(cov.loc[vn1, vn1])

        plt.clf()
        plt.axes([0.15, 0.11, 0.8, 0.8])
        plt.grid(True)
        plt.plot(age, traj, '-', color='black', lw=4)
        plt.fill_between(
            age,
            traj - 2 * se * proj,
            traj + 2 * se * proj,
            color='grey',
            alpha=0.5)
        plt.plot([1, maxage], [0, 0], '-', color='black', lw=3)
        plt.xlabel("Age in years", size=15)
        plt.ylabel("Coefficient for %s%s" % (vn, " growth" if growth else ""), size=15)
        plt.title("Control for height and BMI at time of SBP")
        pdf.savefig()

pdf.close()
