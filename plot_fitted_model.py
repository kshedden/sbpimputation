import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from config import *
import os
import sys
import json

dim = 1
method = sys.argv[1]
qtype = 0

effects = ["age_x", "I(age_x ** 2)", "I(age_x ** 3)", "Female:age_x", "lognummeas",
           "Female:I(age_x ** 2)", "Female:I(age_x ** 3)", "Female:lognummeas",
           "PregMo_Use_cen", "I(PregMo_Use_cen ** 2)", "LactMo_Use_cen",
           "I(LactMo_Use_cen ** 2)", "Female"]

pdf = PdfPages("model_plots.pdf")

for vn in allowed_controls:

    fn = os.path.join(method, "dim_%d" % dim, "%s_%d.txt" % (vn, qtype))
    fid = open(fn)
    pa = {}
    for line in fid:
        for ef in effects:
            if line.startswith(ef + " "):
	            line = line[len(ef):]
	            tok = line.strip().split()
	            pa[ef] = float(tok[0])
    fid.close()

    fid = open("%s_cen.json" % vn)
    cen = json.load(fid)

    # Age plots
    ages = np.linspace(10, 20, 20)
    age_x = (ages - 10) / 10
    male_sbp = pa["age_x"]*age_x + pa["I(age_x ** 2)"]*age_x**2 + pa["I(age_x ** 3)"]*age_x**3
    female_sbp = (male_sbp + pa["Female:age_x"]*age_x + pa["Female:I(age_x ** 2)"]*age_x**2 +
                  pa["Female:I(age_x ** 3)"]*age_x**3 + pa["Female"])

    plt.clf()
    plt.axes([0.13, 0.12, 0.67, 0.8])
    plt.title(vn)
    plt.plot(ages, female_sbp, label="Female", lw=4)
    plt.plot(ages, male_sbp, label="Male", lw=4)
    plt.ylabel("SBP change", size=16)
    plt.xlabel("Age", size=16)
    plt.grid(True)
    ha, lb = plt.gca().get_legend_handles_labels()
    leg = plt.figlegend(ha, lb, "center right")
    leg.draw_frame(False)
    pdf.savefig()

    # lognummeas plots
    nummeas = np.arange(1, 10)
    lognummeas = np.log10(nummeas)
    male_lnm = pa["lognummeas"] * lognummeas
    female_lnm = male_lnm + pa["Female:lognummeas"] * lognummeas
    plt.clf()
    plt.axes([0.13, 0.13, 0.67, 0.8])
    plt.title(vn)
    plt.plot(nummeas, female_lnm, label="Female", lw=4)
    plt.plot(nummeas, male_lnm, label="Male", lw=4)
    plt.ylabel("SBP change", size=16)
    plt.xlabel("Number of SBP measurements", size=16)
    plt.grid(True)
    ha, lb = plt.gca().get_legend_handles_labels()
    leg = plt.figlegend(ha, lb, "center right")
    leg.draw_frame(False)
    pdf.savefig()

    # pregmo plots
    pregmo = np.arange(0, 9)
    pregmo_cen = pregmo - cen["PregMo"]
    pregmo_cen2 = pregmo_cen**2
    female_pregmo = pa["PregMo_Use_cen"] * pregmo_cen + pa["I(PregMo_Use_cen ** 2)"] * pregmo_cen2
    plt.clf()
    plt.title(vn)
    plt.plot(pregmo, female_pregmo, '-', lw=4)
    plt.xlabel("Month of pregnancy", size=16)
    plt.ylabel("SBP change", size=16)
    plt.grid(True)
    pdf.savefig()

    # lactmo plots
    lactmo = np.arange(0, 9)
    lactmo_cen = pregmo - cen["PregMo"]
    lactmo_cen2 = pregmo_cen**2
    female_lactmo = pa["PregMo_Use_cen"] * pregmo_cen + pa["I(PregMo_Use_cen ** 2)"] * pregmo_cen2
    plt.clf()
    plt.title(vn)
    plt.plot(lactmo, female_pregmo, '-', lw=4)
    plt.xlabel("Month of pregnancy", size=16)
    plt.ylabel("SBP change", size=16)
    plt.grid(True)
    pdf.savefig()

pdf.close()
