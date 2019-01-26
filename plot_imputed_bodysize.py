"""
Plot a subset of the imputed trajectories together with the actual data.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from data_tools import get_data
import pandas as pd
import os
from config import *

def plot(impvar):

    # Plot only a few curves to avoid overplotting
    di = []
    for j in range(5):
        di.append(pd.read_csv(os.path.join("imputed_data", "%s_imp_%d.csv" % (impvar, j))))

    dm = []
    for j in range(5):
        pm = di[j].set_index("ID").T.to_dict("series")
        dm.append(pm)

    # Ages of imputed data
    cols = di[0].columns
    cols = [x for x in cols if x.startswith(impvar)]
    ages = np.arange(1, len(cols)+1)

    for female in False, True:

        dx = get_data(female, impvar)

        jj = 0
        for idx, v in dx.groupby("ID"):

            if idx not in dm[0]:
                continue

            plt.clf()

            # Plot the imputed data
            for j in range(5):
                plt.plot(ages, dm[j][idx][0:len(ages)], color='purple',
                         lw=4, alpha=0.5)

            # The mean is the same, only plot it once
            plt.plot(ages, dm[0][idx][len(ages):], color='blue', lw=4)
            plt.plot(v.Age, v[impvar], 'o', color='orange')
            plt.grid(True)
            plt.title("%d" % idx)
            plt.xlabel("Age", size=16)
            plt.ylabel(impvar, size=16)
            pdf.savefig()

            jj += 1
            if jj > 20:
                break

pdf = PdfPages("plot_imputed.pdf")
for iv in allowed_controls:
    plot(iv)
pdf.close()