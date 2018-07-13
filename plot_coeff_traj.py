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

ages = np.arange(1, maxage+1)

for vn in allowed_outcomes:

    f = open("imp_%s.txt" % vn)
    neq = 0
    for line in f:
        if line.startswith("="):
            neq += 1
            if neq == 2:
                break

    next(f)
    rv = []
    for line in f:
        rv.append([float(x) for x in line.rstrip().split()])
    rv = pd.DataFrame(rv, columns=["index", "age", "coeff", "se"])

    aq = np.linspace(-1, 1, len(ages))

    plt.clf()
    plt.grid(True)
    plt.plot(rv.age, rv.coeff, '-', color='orange', lw=4)
    plt.fill_between(rv.age, rv.coeff-2*rv.se, rv.coeff+2*rv.se, color='grey', alpha=0.5)
    plt.xlabel("Age", size=15)
    plt.ylabel("Coefficient for %s" % vn, size=15)
    pdf.savefig()

pdf.close()