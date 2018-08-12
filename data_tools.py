import pandas as pd
import numpy as np
from config import *

# Read the cohort file
pa = "/nfs/kshedden/Beverly_Strassmann/Cohort_2018_080118.csv.gz"
df = pd.read_csv(pa)

# Modify some of the variable names
df["Sex"] = df["Sex"].replace({0: "Female", 1: "Male"})
df["Female"] = df["Sex"].replace({"Female": 1, "Male": 0})
df["Age"] = np.around(df.Age_Yrs, 3)
df["HAZ"] = df.HAZ_18
df["BAZ"] = df.BAZ_18
df["WAZ"] = df.WAZ_18
df["HT"] = df["Ht_Ave_18"]

df["log2T_use"] = np.log(df.T_use) / np.log(2)
df["log2T_use_Z"] = (df.log2T_use - df.log2T_use.mean()) / df.log2T_use.std()

df["Breast_Stage_Z"] = (df.Breast_Stage - df.Breast_Stage.mean()) / df.Breast_Stage.std()

age_mean = df.Age.mean()
age_sd = df.Age.std()
df["AgeZ"] = (df.Age - age_mean) / age_sd

def get_data(female, impvar, others=None):
    """
    Returns the data to be used for imputing 'impvar', either for females
    or for males.  Returns a data frame with ID, Age, AgeZ, and any variables
    listed in 'others'.
    """

    vars = ["ID", "Age", "AgeZ", impvar]
    if others is not None:
        vars.extend(others)

    return df.loc[df.Female == female, vars]

