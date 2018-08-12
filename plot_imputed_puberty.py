"""
Plot a subset of the imputed trajectories together with the actual data.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from data_tools import get_data, df
import pandas as pd
import os

#impvar = "Breast_Stage_Z"
impvar = "log2T_use_Z"

pdf = PdfPages("plot_imputed_%s.pdf" % impvar)

# Plot only a few curves to avoid overplotting
di = []
for j in range(5):
    di.append(pd.read_csv(os.path.join("imputed_data_puberty", "%s_imp_%d.csv" % (impvar, j))))

for female in False, True:

    dx = get_data(female, impvar)

    idx = dx.ID.unique().astype(np.int).tolist()

    jj = 0
    for id0 in idx:

        vv = df.loc[df.ID == id0, :]
        v0 = dx.loc[dx.ID == id0, :]

        plt.clf()
        plt.title("ID=%d" % id0)
        plt.grid(True)

        plt.plot(vv.Age, vv[impvar], 'o', color='purple')

        noplot = False
        for j in range(5):

            v1 = di[j].loc[di[j].ID == id0, :]

            if v1.size == 0:
                noplot = True
                break

            plt.plot(v1.Age, v1[impvar], 'o', color='orange',
                     lw=4, alpha=0.5)

        if noplot:
            continue

        plt.xlabel("Age", size=15)
        plt.ylabel(impvar, size=15)
        pdf.savefig()

        jj += 1
        if jj > 20:
            break

pdf.close()