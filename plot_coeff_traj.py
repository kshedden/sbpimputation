"""
Generate plots of the time-varying coefficients for stature variables
as functions of age.
"""

import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from config import *

pdf = PdfPages("coeff_traj.pdf")

ages = np.arange(1, maxage + 1)

for vn in allowed_controls:
    for jx in 1,2:

        f = open("imp_%s_%d.txt" % (vn, jx))
        cf = []
        for line in f:
            if line.startswith("BEGIN-TRAJ"):
                cx = []
                for line in f:
                    if line.startswith("END-TRAJ"):
                        break
                    cx.append(line.rstrip().split())
                cols = cx[0]
                cx = cx[1:]
                cx = pd.DataFrame(cx, columns=cols)
                for c in cx.columns:
                    cx[c] = pd.to_numeric(cx[c])
                cf.append(cx)

        if "Z" not in vn:
            st = ["change in", "relative change in", "absolute"]
        else:
            # No relative change for Z-score variables
            st = ["change in", "absolute"]

        for jr, rv in enumerate(cf):

            plt.clf()
            plt.grid(True)
            plt.plot(rv.age, rv.coeff, '-', color='orange', lw=4)
            plt.fill_between(
                rv.age,
                rv.coeff - 2 * rv.se,
                rv.coeff + 2 * rv.se,
                color='grey',
                alpha=0.5)
            plt.plot([1, maxage], [0, 0], '-', color='black', lw=4)
            plt.xlabel("Age", size=15)
            plt.ylabel("Coefficient for %s %s" % (st[jr], vn), size=15)
            if jx == 1:
                plt.title("Control for height/BMI at time of SBP")
            else:
                plt.title("Control for %s at time of SBP" % vn)
            pdf.savefig()

pdf.close()
