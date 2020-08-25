"""
Generate estimates of the mean and SE for contrasts comparing
people exposed to different patterns of undernourishment.
"""

import pandas as pd
import numpy as np
import sys, os
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from scipy.stats.distributions import norm
from config import *

bp_var = sys.argv[1]
if len(sys.argv) != 6 or bp_var not in ("sbp", "dbp"):
    print("usage: growth_patterns.py [sbp|dbp] dim pcs_idx [none|zf|zm] [dadbp|nodadbp]\n")
    sys.exit(0)
bp_dir = bp_var

ndim = int(sys.argv[2])

# Growth pattern to compare to reference
pcs_idx = int(sys.argv[3])

# If not "none", change adult height or BMI by one Z-score, not
# by a fixed reference amount.  Set to either "zf" to adjust
# by female Z-scores, or "zm" to adjust by male Z-scores.
zs = sys.argv[4]
if zs not in ("none", "zf", "zm"):
	1/0

dadbp = sys.argv[5] == "dadbp"

sfx = "_1_%d" % pcs_idx

if zs != "none":
    sfx = sfx + "_" + zs

ages = np.arange(1, maxage + 1)

fn = os.path.join(bp_var, "mixed", "dim_%d" % ndim, "growth_patterns_controlcbs_nodadbp%s.txt" % sfx)
if dadbp:
    fn = fn.replace("nodadbp", "dadbp")
fid = open(fn, "w")
fid.write("```\n")

for vn in allowed_controls:

    # This file only handles controlcbs, growth_patterns_nocontrol.py handles the
    # the uncontrolled model.
    fn = os.path.join(bp_dir, "mixed", "dim_%d" % ndim, "%s_controlcbs_nodadbp.txt" % vn)
    if dadbp:
        fn = fn.replace("nodadbp", "dadbp")
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

    # The PC loadings for the childhood exposure variable (based on mean within-imputation
    # covariance)
    cf[0] = [[float(y) for y in x.split(",")] for x in cf[0]]
    proj = np.asarray(cf[0])

    # The (mean within-imputation) variances of the childhood growth
    # variable at each age
    cf[1] = [float(x) for x in cf[1]]
    vax = np.asarray(cf[1])

    from io import StringIO

    params = "\n".join(cf[2])
    params = pd.read_csv(StringIO(params))

    cov = "\n".join(cf[3])
    cov = pd.read_csv(StringIO(cov), index_col=0)

    pa = genpat(len(vax), np.sqrt(vax))

    # HT_cen
    # BMI_cen
    # ?_pc0            ?=childhood exposure variable
    # HT_cen:?_pc0
    # BMI_cen:?_pc0

    for jx in 0, 1:
        if jx == 0:
            # Control for current Ht
            ivn = "Ht" # an adult body dimension
            vu = "cm"  # units for ivn
            if zs == "none":
                hxs = 5
                hxn = str(hxs) + vu
            elif zs == "zf":
                hxs = 5.5
                hxn = "1SD(F)"
            elif zs == "zm":
                hxs = 7.1
                hxn = "1SD(M)"
            else:
                1/0
            v0 = np.r_[1, 0, 0, 0, 0]
            v1 = np.r_[0, 0, 1, 0, 0]
            v2 = np.r_[0, 0, 0, 1, 0]
        elif jx == 1:
            # Control for current BMI
            ivn = "BMI"
            vu = "kg/m^2"
            if zs == "none":
                hxs = 3
                hxn = str(hxs) + vu
            elif zs == "zf":
                hxs = 2.5
                hxn = "1SD(F)"
            elif zs == "zm":
                hxs = 1.8
                hxn = "1SD(M)"
            else:
                1/0
            v0 = np.r_[0, 1, 0, 0, 0]
            v1 = np.r_[0, 0, 1, 0, 0]
            v2 = np.r_[0, 0, 0, 0, 1]
        else:
            1/0

        fid.write(vn + "\n")
        vx, vy, vs = [], [], []

        # The reference point against which comparisons are made
        pcs_ref = np.dot(pa[pcs_idx][1], proj)
        hx_ref = 0
        u0 = hx_ref*v0 + pcs_ref*v1 + hx_ref*pcs_ref*v2

        for hx in -hxs, 0, hxs:
            for x in pa:
                vx.append(x[0])

                # Project the test function against the PC
                pcs = np.dot(x[1], proj)

                u = hx*v0 + pcs*v1 + hx*pcs*v2
                u -= u0

                vy.append(np.dot(u, params.iloc[:, 1]))
                vs.append(np.dot(u, np.dot(cov, u)))

        vy = np.reshape(np.asarray(vy), (-1, 7)).T
        vs = np.reshape(np.asarray(vs), (-1, 7)).T
        vs = np.sqrt(vs)
        col = ["%s:-%s" % (ivn, hxn), "%s:0" % ivn, "%s:+%s" % (ivn, hxn)]
        vy = pd.DataFrame(vy, index=vx[0:7], columns=col)
        fid.write(vy.to_string())
        fid.write("\n\n" + vn + " Standard errors:\n")
        vs = pd.DataFrame(vs, index=vx[0:7], columns=col)
        fid.write(vs.to_string())
        vz = vy / vs
        vp = 2 * norm.cdf(-np.abs(vz))
        vp = pd.DataFrame(vp, index=vz.index, columns=vz.columns)
        fid.write("\n\n" + vn + " p-values:\n")
        fid.write(vp.to_string())
        fid.write("\n\n")

fid.write("```\n")
fid.close()
