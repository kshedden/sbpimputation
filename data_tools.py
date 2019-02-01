import pandas as pd
import numpy as np
from config import *

log = open("data_tools.log", "w")

# Read the cohort file
pa = "/nfs/kshedden/Beverly_Strassmann/Cohort_2018.csv.gz"
df = pd.read_csv(pa)

log.write("%d distinct subjects\n" % df.ID.unique().size)
df = df.loc[df.Bamako.isin([0, 1]), :]
log.write("%d after requiring Bamako in (0, 1)\n" % df.ID.unique().size)
df = df.loc[df.Ht_Lag_Flag!=4, :]
log.write("%d after requiring Ht_Lag_Flag!=4\n" % df.ID.unique().size)
df = df.loc[df.survival==1, :]
log.write("%d after requiring survival=1\n" % df.ID.unique().size)

# Modify some of the variable names
df["Sex"] = df["Sex"].replace({0: "Female", 1: "Male"})
df["Female"] = df["Sex"].replace({"Female": 1, "Male": 0})
df["Age"] = np.around(df.Age_Yrs, 3)
df["HAZ"] = df.HAZ_18
df["BAZ"] = df.BAZ_18
df["WAZ"] = df.WAZ_18
df["BMI"] = df.BMI_18
df["HT"] = df["Ht_Ave_18"]

df["log2T_use"] = np.log(df.T_use) / np.log(2)
df["log2T_use_Z"] = (df.log2T_use - df.log2T_use.mean()) / df.log2T_use.std()

df["Breast_Stage_Use_Z"] = (df.Breast_Stage_Use - df.Breast_Stage_Use.mean()) / df.Breast_Stage_Use.std()

age_mean = df.Age.mean()
age_sd = df.Age.std()
df["AgeZ"] = (df.Age - age_mean) / age_sd

df["Pregnant"] = df.ReproStat.isin([2, 4])
df["Lactating"] = df.ReproStat.isin([3])

# Require people to have at least one SBP measurement
idx = df.loc[pd.notnull(df.SBP_MEAN), :].groupby("ID").head(1)["ID"]
df = df.loc[df.ID.isin(idx), :]
log.write("%d after requiring at least one SBP measurement\n" % df.ID.unique().size)

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
