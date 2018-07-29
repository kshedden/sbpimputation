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

for vn in allowed_outcomes:

    f = open("imp_%s.txt" % vn)
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

    for deriv in 0, 1:

        rv = cf[deriv]

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
        plt.ylabel("Coefficient for %s %s" % (["change in", ""][deriv], vn), size=15)
        pdf.savefig()

pdf.close()
