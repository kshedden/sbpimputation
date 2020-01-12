import pandas as pd
from data_tools import df
import os
import numpy as np

dq = df.loc[:, ["ID", "Female", "Age", "BMI_18", "WT", "HT"]]

dq["Age"] = np.around(dq.Age)

def fl(x):
    return x.dropna().size

su = dq.groupby(["Age", "Female"])["BMI_18"].agg([fl, np.mean, np.std])
su = su.rename(columns={"fl": "n"})
su["n"] = su["n"].astype(np.int)
su.to_csv("bmi_stats.csv")

su = dq.groupby(["Age", "Female"])["WT"].agg([fl, np.mean, np.std])
su = su.rename(columns={"fl": "n"})
su["n"] = su["n"].astype(np.int)
su.to_csv("wt_stats.csv")

su = dq.groupby(["Age", "Female"])["HT"].agg([fl, np.mean, np.std])
su = su.rename(columns={"fl": "n"})
su["n"] = su["n"].astype(np.int)
su.to_csv("ht_stats.csv")
