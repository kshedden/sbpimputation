"""
Plot a subset of the imputed trajectories together with the actual data.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from data_tools import get_data
import pandas as pd
import os

impvar = "BAZ"

pdf = PdfPages("plot_imputed_%s.pdf" % impvar)

# Plot only a few curves to avoid overplotting
di = []
for j in range(5):
    di.append(pd.read_csv(os.path.join("imputed_data", "%s_imp_%d.csv" % (impvar, j))))

dm = []
for j in range(5):
    pm = {}
    for i in range(di[j].shape[0]):
        row = di[j].iloc[i, :]
        pm[row.ID] = np.asarray(row.iloc[1:])
    dm.append(pm)

# Get the ages to plot
ages = []
for x in di[0].columns:
    if x.startswith(impvar):
        a = int(x[len(impvar):])
        ages.append(a)
ages = np.asarray(ages)

for female in False, True:

    dx = get_data(female, impvar)

    jj = 0
    for idx, v in dx.groupby("ID"):
        plt.clf()

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

pdf.close()