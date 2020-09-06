import sys
import os

# Run models for all outcomes in parallel
#
# python run_sbp_pca.py [sbp|dbp] [ndim] [controlcbs|nocontrolcbs] [dadbp|nodadbp]

# Replaces sbp_pca.sh

vn = ["HT", "HAZ", "BAZ", "WAZ", "WT", "BMI"]

if len(sys.argv) != 6:
    print("Wrong arguments")
    1/0

bptype = sys.argv[1]
if bptype not in ["sbp", "dbp"]:
    raise ValueError("bptype must be sbp or dbp")

ndim = int(sys.argv[2])

cntrl_type = sys.argv[3]
if cntrl_type not in ["controlcbs", "nocontrolcbs"]:
    raise ValueError("Invalid control type")

dadbp = sys.argv[4]
if dadbp not in ["dadbp", "nodadbp"]:
    raise ValueError("Invalid dadbp")

c = []
for v in vn:
    for sex in "female", "male":
        s = "%s %s %d %s %s %s" % (v, bptype, ndim, cntrl_type, dadbp, sex)
        c.append(s)

cmd = ":".join(c)

cmd = "echo %s | rush -D \":\" -k 'python sbp_pca_sex.py {}'" % cmd

os.system(cmd)







