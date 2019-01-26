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
    print("usage: growth_patterns.py [sbp|dbp] dim\n")
    sys.exit(0)
bp_dir = bp_var

ndim = int(sys.argv[2])

sfx = ""

# Could make this more general
fid = open("%s/mixed/dim_%d/growth_patterns%s.txt" % (bp_var, ndim, sfx))
pdf = PdfPages("%s/mixed/dim_%d/growth_patterns_plot%s.pdf" % (bp_var, ndim, sfx))

lb = ["1→1", "0→0", "-1→0", "0→-1", "-1→-1"]

cols = ["#003f5c", "#bc5090", "#ff6361", "#58508d", "#ffa600"]

next(fid) # Skip ```

fq = norm.ppf(1 - 0.025/14)

while True:

    group = next(fid).rstrip()
    if group == "```":
        break
    head = next(fid).strip().split()
    table = [next(fid).rstrip() for k in range(5)]
    rows = [x[0:35].rstrip() for x in table]
    table = [x[35:].lstrip() for x in table]
    next(fid)
    next(fid)
    next(fid)
    se_table = [next(fid).rstrip() for k in range(5)]
    se_table = [x[35:].lstrip() for x in se_table]
    next(fid)

    head = [re.sub("(:)([1-9])", r":+\2", x) for x in head]

    table = pd.read_csv(StringIO("\n".join(table)), delim_whitespace=True, header=None)
    table.columns = head
    table.index = rows

    se_table = pd.read_csv(StringIO("\n".join(se_table)), delim_whitespace=True, header=None)
    se_table.columns = head
    se_table.index = rows

    if group not in ("HT", "BMI"):
        print("Skipping %s" % group)
        continue
    print("Processing %s" % group)

    table = table.unstack()
    se_table = se_table.unstack()
    tab = pd.concat((table, se_table), axis=1)
    tab.columns = ["Est", "SE"]

    plt.figure(figsize=(10, 5))
    plt.clf()
    ax = plt.axes([0.1, 0.15, 0.88, 0.75])
    plt.plot([-1, 15], [0, 0], '-', color='grey')
    jj = np.arange(tab.shape[0])
    for j in jj:
        plt.plot([j, j], [tab.Est[j]-fq*tab.SE[j], tab.Est[j]+fq*tab.SE[j]], color=cols[j%5])
        if tab.SE[j] != 0:
            plt.plot([j], [tab.Est[j]], 'o', color=cols[j%5], ms=5)
        else:
            plt.plot([j], [tab.Est[j]], 's', mfc="none", mec=cols[j%5], ms=5)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    plt.ylabel("Change in %s (mm/Hg)" % bp_var.upper(), size=16)
    plt.figtext(0.29, 0.04, head[0], ha='center')
    plt.figtext(0.54, 0.04, head[1], ha='center')
    plt.figtext(0.78, 0.04, head[2], ha='center')
    plt.title("Childhood %s trajectory" % group)

    ax.set_xticks(jj)
    ax.set_xticklabels([x for x in lb + lb + lb])

    for j in jj:
        ax.get_xticklabels()[j].set_color(cols[j%5])

    plt.xlim(-0.5, 14.5)

    pdf.savefig()

pdf.close()

