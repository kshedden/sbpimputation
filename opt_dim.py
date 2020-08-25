"""
Prepare a table of AIC values for models fit to different dimensional
representiations of the early childhood growth trajectory.
"""

import sys, os
import numpy as np

bp_var = sys.argv[1]
if bp_var not in ("sbp", "dbp"):
    print("usage: python opt_dim.py [sbp|dbp]\n")
    sys.exit(0)
bp_dir = bp_var

vars = ["HT", "WT", "BMI", "HAZ", "WAZ", "BAZ"]
dims = [1, 2]

out = open(os.path.join(bp_dir, "mixed", "opt_dim.txt"), "w")
out.write("```\n")
out.write("          Outcome                            LL                          AIC                   Dimension\n")
for va in vars:
    for controlcbs in False, True:
        for dadbp in False, True:

            ll = []
            for d in dims:

                fn = os.path.join(bp_dir, "mixed", "dim_%d" % d, "%s_nocontrolcbs_nodadbp.txt" % va)
                if controlcbs:
                    fn = fn.replace("nocontrolcbs", "controlcbs")
                if dadbp:
                    fn = fn.replace("nodadbp", "dadbp")
                fid = open(fn)
                for line in fid:
                    if line.startswith("mean IC"):
                        ll.append(float(line.split()[2]))
                        break

            if len(ll) != len(dims):
                1/0

            ll = np.asarray(ll)

            # Convert to AIC
            aic = ll - 3 * np.arange(len(dims))

            icm = np.argmax(aic)

            l1 = "control" if controlcbs else "no control"
            l2 = "dadbp" if dadbp else "no dadbp"

            out.write("%-5s %-12s %-10s %-20s %-20s %10d\n" % (va, l1, l2, ll, aic, icm + 1))

out.write("```\n")
out.close()
