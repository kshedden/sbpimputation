import sys
import os

# Run models for all outcomes in parallel
#
# python run_sbp_pca.py [sbp|dbp] [ndim]

# Replaces sbp_pca.sh

#vn = ["HT", "HAZ", "BAZ", "WAZ", "WT", "BMI"]
vn = ["HT", "WT"]

if len(sys.argv) != 3:
    print("Wrong arguments")

bptype = sys.argv[1]
if bptype not in ["sbp", "dbp"]:
    raise ValueError("bptype must be sbp or dbp")

ndim = int(sys.argv[2])

c = []
for v in vn:
    for sex in "female", "male":
        for adj in 0, 1, 2, 3, 4:
            s = "%s %s %d %s %d" % (v, bptype, ndim, sex, adj)
            c.append(s)

cmd = ":".join(c)

cmd = "echo %s | rush -D \":\" -k 'python sbp_pca_sex.py {}'" % cmd

os.system(cmd)
