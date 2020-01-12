"""
Generate a plot showing mean BP differences between people
exposed to different patterns of undernourishment.
"""

import sys
sys.path.insert(0, "/afs/umich.edu/user/k/s/kshedden/statsmodels_fork/statsmodels")
import numpy as np
import matplotlib.pyplot as plt
from io import StringIO
import pandas as pd
import re
from scipy.stats.distributions import norm
from matplotlib.backends.backend_pdf import PdfPages

bp_var = sys.argv[1]
if bp_var not in ("sbp", "dbp"):
    print("usage: growth_patterns_plot.py [sbp|dbp] dim\n")
    sys.exit(0)
bp_dir = bp_var

ndim = int(sys.argv[2])

# Select the reference category
sfx = "_1_3"

# Could make this more general
fname = "%s/mixed/dim_%d/growth_patterns%s.txt" % (bp_var, ndim, sfx)
print("Reading %s\n" % fname)
fid = open(fname)

pdf = PdfPages("%s/mixed/dim_%d/growth_patterns_plot%s.pdf" % (bp_var, ndim, sfx))

lb = ["1→1", "0→1", "1→0", "0→0", "-1→0", "0→-1", "-1→-1"]

from turbo_colormap import turbo_colormap_data
c = turbo_colormap_data
ii = np.linspace(15, len(c) - 30, 7) # Avoid very dark colors
ii = np.round(ii).astype(np.int)
cols = [c[i] for i in ii]
syms = ['s', 'o', 'x', '+', 'D', '>', '<']

next(fid) # Skip initial ```

fq = norm.ppf(1 - 0.025/14)

while True:

    group = next(fid).rstrip()
    if group == "```":
        break

    head = next(fid).strip().split()
    table = [next(fid).rstrip() for k in range(7)]

    rows = [x[0:36].rstrip() for x in table]
    table = [x[36:].lstrip() for x in table]
    next(fid)
    next(fid)
    next(fid)
    se_table = [next(fid).rstrip() for k in range(7)]
    se_table = [x[36:].lstrip() for x in se_table]
    next(fid)
    next(fid)
    next(fid)
    pv_table = [next(fid).rstrip() for k in range(7)]
    pv_table = [x[36:].lstrip() for x in se_table]
    next(fid)

    agrp = head[0].split(":")[0]
    head = [x.split(":")[1] for x in head]
    head = [re.sub("^([1-9])", r"+\1", x) for x in head]

    table = pd.read_csv(StringIO("\n".join(table)), delim_whitespace=True, header=None)
    table.columns = head
    table.index = rows

    se_table = pd.read_csv(StringIO("\n".join(se_table)), delim_whitespace=True, header=None)
    se_table.columns = head
    se_table.index = rows

    if group not in ("HT", "BMI", "WT"):
        print("Skipping %s" % group)
        continue
    print("Processing %s" % group)

    table = table.unstack()
    se_table = se_table.unstack()
    tab = pd.concat((table, se_table), axis=1)
    tab.columns = ["Est", "SE"]

    plt.figure(figsize=(6, 3.5))
    plt.clf()
    ax = plt.axes([0.18, 0.17, 0.8, 0.7])
    plt.plot([-1, 25], [0, 0], '-', color='grey')
    jj = np.r_[0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22]
    for i, j in enumerate(jj):
        if j < len(syms):
            plt.plot([j, j], [tab.Est[i]-fq*tab.SE[i], tab.Est[i]+fq*tab.SE[i]],
                color=cols[i%len(syms)], lw=4, label=lb[j])
        else:
            plt.plot([j, j], [tab.Est[i]-fq*tab.SE[i], tab.Est[i]+fq*tab.SE[i]],
                color=cols[i%len(syms)], lw=4)
#        if tab.SE[i] != 0:
#            dl = {"mfc": "none", "color": cols[i%len(syms)], "ms": 7}
#            if j < len(syms):
#                dl["label"] = lb[j]
#            plt.plot([j], [tab.Est[i]], syms[i%len(syms)], **dl)
#        else:
#            plt.plot([j], [tab.Est[i]], 's', mfc="none", mec=cols[i%len(syms)], ms=5)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    plt.ylabel("%s difference (mm Hg)" % bp_var.upper(), size=12)
    plt.figtext(0.01, 0.01, "Adult %s\ndifference:" % agrp.replace("Ht", "height"),
                ha="left", va='bottom', size=12)
    plt.figtext(0.01, 0.93, "Childhood %s z:" % group.replace("HT", "height"), ha="left", size=12)

    ax.set_xticks([3, 11, 19])
    for j in range(len(head)):
        if "kg" in head[j]:
            head[j] = head[j].replace("kg/m^2", "{\sf kg}/{\sf m}^2")
            head[j] = "$" + head[j] + "$"
    ax.set_xticklabels(head)

    plt.xlim(-0.5, 22.5)

    ha, lbx = plt.gca().get_legend_handles_labels()
    leg = plt.figlegend(ha, lbx, (0.3, 0.87), ncol=4)
    leg.draw_frame(False)

    for x in 7, 15:
        plt.plot([x, x], [-10, 10], '-', color='grey')

    plt.plot([11,], [0,], '|', ms=5, color='black')

    plt.ylim(-10, 10)
    pdf.savefig()

pdf.close()
