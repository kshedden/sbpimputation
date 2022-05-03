import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from config import *
import os
import sys
import json

bp_var = sys.argv[1]
if bp_var not in ("sbp", "dbp"):
    print("usage: plot_fitted_model.py [sbp|dbp] [mixed|gee]\n")
    sys.exit(0)
bp_dir = bp_var

dim = 1
tr = "mixed"

effects = ["age_x", "I(age_x ** 2)", "I(age_x ** 3)", "Female:age_x", "lognummeas",
           "Female:I(age_x ** 2)", "Female:I(age_x ** 3)", "Female:lognummeas",
           "PregMo_Use_cen", "I(PregMo_Use_cen ** 2)", "LactMo_Use_cen",
           "I(LactMo_Use_cen ** 2)", "Female", "HT_cen", "BMI_cen"]

pdf = PdfPages(os.path.join(bp_dir, tr, "dim_%d" % dim, "model_plots.pdf"))

for vn in allowed_controls:
    for controlcbs in False, True:
        for dadbp in False, True:

            fn = os.path.join(bp_dir, tr, "dim_%d" % dim, "%s_nocontrolcbs_nodadbp_adj_female.txt" % vn)
            if controlcbs:
                fn = fn.replace("nocontrolcbs", "controlcbs")
            if dadbp:
                fn = fn.replace("nodadbp", "dadbp")
            fid = open(fn)
            pa = {}
            for line in fid:
                for ef in effects:
                    if line.startswith(ef + " "):
        	            line = line[len(ef):]
        	            tok = line.strip().split()
        	            pa[ef] = float(tok[0])
            fid.close()

            fid = open(os.path.join(bp_dir, tr, "dim_%d" % dim, "%s_cen.json" % vn))
            cen = json.load(fid)

            # Age plots
            ages = np.linspace(10, 20, 20)
            age_x = (ages - 10) / 10
            male_sbp = pa["age_x"]*age_x + pa["I(age_x ** 2)"]*age_x**2 + pa["I(age_x ** 3)"]*age_x**3
            female_sbp = (male_sbp + pa["Female:age_x"]*age_x + pa["Female:I(age_x ** 2)"]*age_x**2 +
                          pa["Female:I(age_x ** 3)"]*age_x**3 + pa["Female"])
            mm = (np.mean(male_sbp) + np.mean(female_sbp)) / 2
            male_sbp -= mm
            female_sbp -= mm
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
            mm = (np.mean(male_lnm) + np.mean(female_lnm)) / 2
            male_lnm -= mm
            female_lnm -= mm
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
            female_pregmo -= female_pregmo.mean()
            plt.clf()
            plt.title(vn)
            plt.plot(pregmo, female_pregmo, '-', lw=4)
            plt.xlabel("Month of pregnancy", size=16)
            plt.ylabel("SBP change", size=16)
            plt.grid(True)
            pdf.savefig()

            # lactmo plots
            lactmo = np.arange(0, 24)
            lactmo_cen = lactmo - cen["LactMo"]
            lactmo_cen2 = lactmo_cen**2
            female_lactmo = pa["LactMo_Use_cen"] * lactmo_cen
            female_lactmo -= female_lactmo.mean()
            plt.clf()
            plt.title(vn)
            plt.plot(lactmo, female_lactmo, '-', lw=4)
            plt.xlabel("Month of lactation", size=16)
            plt.ylabel("SBP change", size=16)
            plt.grid(True)
            pdf.savefig()

            # Current BMI plots
            if controlcbs:
                bmi = np.arange(12, 25)
                bmi_cen = bmi - cen["BMI"]
                bmi_effect = pa["BMI_cen"] * bmi_cen
                bmi_effect -= bmi_effect.mean()
                plt.clf()
                plt.axes([0.13, 0.13, 0.67, 0.8])
                plt.title(vn)
                plt.plot(bmi, bmi_effect, '-', lw=4)
                plt.xlabel("Adult BMI", size=16)
                plt.ylabel("SBP change", size=16)
                plt.grid(True)
                leg.draw_frame(False)
                pdf.savefig()

            # Current height plots
            if controlcbs:
                hgt = np.arange(120, 183)
                hgt_cen = hgt - cen["HT"]
                hgt_effect = pa["HT_cen"] * hgt_cen
                hgt_effect -= hgt_effect.mean()
                plt.clf()
                plt.axes([0.13, 0.13, 0.67, 0.8])
                plt.title(vn)
                plt.plot(hgt, hgt_effect, '-', lw=4)
                plt.xlabel("Adult height", size=16)
                plt.ylabel("SBP change", size=16)
                plt.grid(True)
                pdf.savefig()

pdf.close()
