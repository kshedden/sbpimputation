import pandas as pd

maxage = 11

# Read the cohort file
pa = "/nfs/kshedden/Beverly_Strassmann/Cohort_2018.csv.gz"
df = pd.read_csv(pa)

# Modify some of the variable names
df["Sex"] = df["Sex"].replace({0: "Female", 1: "Male"})
df["Female"] = df["Sex"].replace({"Female": 1, "Male": 0})
df["Age"] = df.Age_Yrs
df["HAZ"] = df.HAZ_17
df["BAZ"] = df.BAZ_15
df["WAZ"] = df.WAZ_15
df["HT"] = df["Ht_Ave"]

df = df.loc[df.Age <= maxage, :]


def get_data(female, impvar):
    """
    Returns the data to be used for imputing 'impvar', either for females
    or for males.
    """
    return df.loc[df.Female == female, ["ID", "Age", impvar]].dropna()
