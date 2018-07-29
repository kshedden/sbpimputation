import pandas as pd
from config import *

# Read the cohort file
pa = "/nfs/kshedden/Beverly_Strassmann/Cohort_2018_072818.csv.gz"
df = pd.read_csv(pa)

# Modify some of the variable names
df["Sex"] = df["Sex"].replace({0: "Female", 1: "Male"})
df["Female"] = df["Sex"].replace({"Female": 1, "Male": 0})
df["Age"] = df.Age_Yrs
df["HAZ"] = df.HAZ_18
df["BAZ"] = df.BAZ_18
df["WAZ"] = df.WAZ_18
df["HT"] = df["Ht_Ave_18"]

# Use data one year beyond the imputation range.
df = df.loc[df.Age <= maxage + 1, :]


def get_data(female, impvar):
    """
    Returns the data to be used for imputing 'impvar', either for females
    or for males.
    """
    return df.loc[df.Female == female, ["ID", "Age", impvar]].dropna()
