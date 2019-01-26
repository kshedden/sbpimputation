"""
Generate estimates of the mean and SE for contrasts comparing
people exposed to different patterns of undernourishment.
"""

import pandas as pd
import numpy as np
import sys, os
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from config import *

bp_var = sys.argv[1]
if bp_var not in ("sbp", "dbp"):
    print("usage: growth_patterns.py [sbp|dbp] dim\n")
    sys.exit(0)
bp_dir = bp_var

ndim = int(sys.argv[2])

# 0, 1, or 2 for reference adult body size being -/0/+ relative to the set point
# Default is 1
hx_idx = 1

# Index of reference child growth exposure, default is 1
pcs_idx = 1

sfx = ""
if hx_idx != 1 or pcs_idx != 1:
    sfx = "_%d_%d" % (hx_idx, pcs_idx)

ages = np.arange(1, maxage + 1)

def genpat(d, sd):
    """
    Generate the test functions.
    """

    pat = []

    x = np.ones(d) * sd
    pat.append(["Always above mean (Z=1 throughout)", x])

    x = np.zeros(d) * sd
    pat.append(["Always at mean (Z=0 throughout)", x])

    x = np.linspace(-1, 0, d) * sd
    pat.append(["Catch up (Z goes from -1 to 0)", x])

    x = np.linspace(0, -1, d) * sd
    pat.append(["Catch down (Z goes from 0 to -1)", x])

    x = -np.ones(d) * sd
    pat.append(["Always below mean (Z=-1 throughout)", x])

    return pat


fid = open(os.path.join(bp_var, "mixed", "dim_%d" % ndim, "growth_patterns%s.txt" % sfx), "w")
fid.write("```\n")

for vn in allowed_controls:
    for growth in False, True:

        if growth and "Z" in vn:
            continue

        fn = os.path.join(bp_dir, "mixed", "dim_%d" % ndim, "%s%s.txt" % (vn, "_growth" if growth else ""))
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
        cf[0] = [float(x) for x in cf[0]]
        proj = np.asarray(cf[0])

        # The (mean within-imputation) variances of the childhood growth
        # variable at each age
        cf[1] = [float(x) for x in cf[1]]
        vax = np.asarray(cf[1])

        from io import StringIO

        params = "\n".join(cf[2])
        params = pd.read_csv(StringIO(params), header=None)

        cov = "\n".join(cf[3])
        cov = pd.read_csv(StringIO(cov), index_col=0)

        pa = genpat(len(vax), np.sqrt(vax))

        # HT_cen
        # BMI_cen
        # Female:BMI_cen
        # ?_pc0            ?=childhood exposure variable
        # HT_cen:?_pc0
        # BMI_cen:?_pc0

        for jx in 0, 1, 2:
            if jx == 0:
                # Ht
                ivn = "Ht" # an adult body dimension
                vu = "cm"  # units for ivn
                hxs = 5    # step size for varying ivn
                v0 = np.r_[1, 0, 0, 0, 0, 0]
                v1 = np.r_[0, 0, 0, 1, 0, 0]
                v2 = np.r_[0, 0, 0, 0, 1, 0]
            elif jx == 1:
                # BMI for males
                ivn = "MBMI"
                vu = ""
                hxs = 3
                v0 = np.r_[0, 1, 0, 0, 0, 0]
                v1 = np.r_[0, 0, 0, 1, 0, 0]
                v2 = np.r_[0, 0, 0, 0, 0, 1]
            elif jx == 2:
                # BMI for females
                ivn = "FBMI"
                vu = ""
                hxs = 3
                v0 = np.r_[0, 1, 1, 0, 0, 0]
                v1 = np.r_[0, 0, 0, 1, 0, 0]
                v2 = np.r_[0, 0, 0, 0, 0, 1]
            else:
                1/0

            vnt = vn + {True: " growth", False: ""}[growth]
            fid.write(vnt + "\n")
            vx, vy, vs = [], [], []

            # The reference point against which comparisons are made
            pcs_ref = np.dot(pa[pcs_idx][1], proj)
            hx_ref = [-hxs, 0, hxs][hx_idx]
            u0 = hx_ref*v0 + pcs_ref*v1 + hx_ref*pcs_ref*v2

            for hx in -hxs, 0, hxs:
                for x in pa:
                    vx.append(x[0])

                    # Project the test function against the PC
                    pcs = np.dot(x[1], proj)

                    u = hx*v0 + pcs*v1 + hx*pcs*v2
                    u -= u0
                    vy.append(np.dot(u, params[1]))
                    vs.append(np.dot(u, np.dot(cov, u)))
            vy = np.reshape(np.asarray(vy), (-1, 5)).T
            vs = np.reshape(np.asarray(vs), (-1, 5)).T
            col = ["%s:-%d%s" % (ivn, hxs, vu), "%s:0%s" % (ivn, vu), "%s:%d%s" % (ivn, hxs, vu)]
            vy = pd.DataFrame(vy, index=vx[0:5], columns=col)
            fid.write(vy.to_string())
            fid.write("\n\n" + vnt + " Standard errors:\n")
            vs = pd.DataFrame(vs, index=vx[0:5], columns=col)
            fid.write(vs.to_string())
            fid.write("\n\n")

fid.write("```\n")
fid.close()
