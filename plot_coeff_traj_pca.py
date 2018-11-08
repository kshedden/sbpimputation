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

bp_dir = "sbp" if sys.argv[1] == "sbp" else "dbp"

mixed = sys.argv[2] == "mixed"
di = "mixed" if mixed else "gee"

ndim = int(sys.argv[3])

pdf = PdfPages(os.path.join(bp_dir, di, "dim_%d" % ndim, "coeff_traj_pca.pdf"))

ages = np.arange(1, maxage + 1)

for vn in allowed_controls:
    for growth in False, True:

        vn_full = vn + ("_growth" if growth else "")

        for varv in 0, 1:

            if growth and "Z" in vn:
                continue

            fn = os.path.join(bp_dir, di, "dim_%d" % ndim, "%s_%d.txt" % (vn_full, varv))
            f = open(fn)

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

            for jr, rv in enumerate(cf):

                if growth:
                    rv.age += 1.5

                plt.clf()
                plt.axes([0.15, 0.11, 0.8, 0.8])
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
                plt.ylabel("Coefficient for %s%s" % (vn, " growth" if growth else ""), size=15)
                if varv == 1:
                    plt.title("Control for height/BMI at time of SBP")
                pdf.savefig()

pdf.close()
