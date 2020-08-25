import pandas as pd
import numpy as np
from config import *
import json

log = open("data_tools.log", "w")

centering = {}

# Read the cohort file
pa = "/nfs/kshedden/Beverly_Strassmann/Cohort_2020.csv.gz"
df = pd.read_csv(pa)
log.write("Original file:\n    %d distinct subjects, %d records\n" %
          (df.ID.unique().size, df.shape[0]))

# The current version of the cohort file has two discordant Sex values,
# fix them here.
def f(x):
    y = x.value_counts()
    return y.index[0]
df["Sex"] = df.groupby("ID")["Sex"].transform(f)

# Modify some of the variable names
df["Sex"] = df["Sex"].replace({0: "Female", 1: "Male"})
df["Female"] = df["Sex"].replace({"Female": 1, "Male": 0})
df["Age"] = np.around(df.Age_Yrs, 3)
df["HAZ"] = df.HAZ
df["BAZ"] = df.BAZ
df["WAZ"] = df.WAZ
df["BMI"] = df.BMI
df["HT"] = df["Ht_Ave_Use"]
df["School"] = df["School_Use"]
df["lognummeas"] = df["lognummeas_BASE10"]
df["Mom_BMI"] = df["F0_Mom_BMI"]
df["Dad_BMI"] = df["F0_Dad_BMI"]
df["MomIdUnique"] = df["F0_MomID.Unique"]
df["Mom_Ht_Ave"] = df["F0_Mom_Ht_Ave"]
df["Dad_Ht_Ave"] = df["F0_Dad_Ht_Ave"]
df["Mom_SBP_Z_Res_USE"] = df["F0_Mom_SBP_Z_Res_USE"]
df["Dad_SBP_Z_Res_USE"] = df["F0_Dad_SBP_Z_Res_USE"]

# Drop people with no Ht data
x = df[["ID", "HT", "Age"]].dropna()
x = x.loc[x.Age <= maxage+1, :]
df = df.loc[df.ID.isin(x.ID), :]
log.write("After dropping people with no HT before age %d:\n" % (maxage + 1))
log.write("    %d distinct subjects, %d records\n" % (df.ID.unique().size, df.shape[0]))

# Drop people with no SBP data
x = df[["ID", "SBP_MEAN"]].dropna()
df = df.loc[df.ID.isin(x.ID), :]
log.write("After dropping people with no SBP\n")
log.write("    %d distinct subjects, %d records\n" % (df.ID.unique().size, df.shape[0]))
df = df.loc[df.Ht_Lag_Flag!=4, :]
log.write("After requiring Ht_Lag_Flag != 4:\n")
log.write("    %d distinct subjects, %d records\n" % (df.ID.unique().size, df.shape[0]))

df = df.loc[df.Bamako.isin([0, 1]), :]
log.write("After requiring Bamako in 0, 1:\n")
log.write("    %d distinct subjects, %d records\n" % (df.ID.unique().size, df.shape[0]))

df["log2T_use"] = np.log(df.T_use) / np.log(2)
m = [df.log2T_use.mean(), df.log2T_use.std()]
df["log2T_use_Z"] = (df.log2T_use - m[0]) / m[1]
centering["log2T_use"] = m

m = [df.Breast_Stage_Use.mean(), df.Breast_Stage_Use.std()]
df["Breast_Stage_Use_Z"] = (df.Breast_Stage_Use - m[0]) / m[1]
centering["Breast_Stage_Use"] = m

age_mean = df.Age.mean()
age_sd = df.Age.std()
df["AgeZ"] = (df.Age - age_mean) / age_sd

df["Pregnant"] = df.ReproStat.isin([2, 4])
df["Lactating"] = df.ReproStat.isin([3])

# Require people to have at least one SBP measurement
idx = df.loc[pd.notnull(df.SBP_MEAN), :].groupby("ID").head(1)["ID"]
df = df.loc[df.ID.isin(idx), :]
log.write("After requiring at least one blood pressure measurement:\n")
log.write("    %d distinct subjects, %d records\n" % (df.ID.unique().size, df.shape[0]))

log.close()

def get_data(female, impvar, others=None):
    """
    Returns the data to be used for imputing 'impvar', either for females
    or for males.  Returns a data frame with ID, Age, AgeZ, and any variables
    listed in 'others'.
    """

    vars = ["ID", "Age", "AgeZ", impvar]
    if others is not None:
        vars.extend(others)

    dx = df.loc[df.Female == female, vars]
    dx["ID"] = dx["ID"].astype(np.int)

    return dx

f = open("centering.json", "w")
json.dump(centering, f)
f.close()
