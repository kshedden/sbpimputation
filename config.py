import numpy as np

# Maximum age for the early childhood window
maxage = 10

# Number of imputated data sets
n_imp = 20

allowed_controls = ("HT", "WT", "HAZ", "WAZ", "BAZ", "BMI")

def genpat(d, sd):
    """
    Generate the test functions.
    """

    pat = []

    x = np.ones(d) * sd
    pat.append(["Always above mean (Z=1 throughout)", x])

    x = np.linspace(0, 1, d) * sd
    pat.append(["High catch up (Z goes from 0 to 1)", x])

    x = np.linspace(1, 0, d) * sd
    pat.append(["High catch down (Z goes from 1 to 0)", x])

    x = np.zeros(d) * sd
    pat.append(["Always at mean (Z=0 throughout)", x])

    x = np.linspace(-1, 0, d) * sd
    pat.append(["Catch up (Z goes from -1 to 0)", x])

    x = np.linspace(0, -1, d) * sd
    pat.append(["Catch down (Z goes from 0 to -1)", x])

    x = -np.ones(d) * sd
    pat.append(["Always below mean (Z=-1 throughout)", x])

    return pat
