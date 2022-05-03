import numpy as np
import pandas as pd
import os

# A class for producing imputed data sets
class mimi(object):

    ix = 0

    def __init__(self, impvar, vx, vb, mn, proj, varls, df, log, bp_var, bp_dir):
        self.impvar = impvar
        self.vx = vx
        self.vb = vb
        self.mn = mn
        self.proj = proj
        self.varls = varls
        self.df = df
        self.log = log
        self.bp_var = bp_var
        self.bp_dir = bp_dir

    def update(self):

        # Create an array to hold the gridded values of impvar
        fn = os.path.join("imputed_data", "%s_imp_%d.csv" % (self.impvar, self.ix))
        di = pd.read_csv(fn)
        di["ID"] = di["ID"].astype(np.int)
        dd = di[self.vx]

        # Get the PC scores for the imputed variable
        dbs = np.dot(dd - self.mn, self.proj)
        dbs = np.concatenate((di.ID.values[:, None], dbs), axis=1)
        dbs = pd.DataFrame(dbs, columns=["ID"] + self.vb)
        dbs.ID = dbs.ID.astype(np.int)

        # Merge all variables created above with other cohort file variables
        dx = pd.merge(
            self.df.loc[:, self.varls], dbs, left_on="ID", right_on="ID", how="left")

        if self.ix == 0:
            self.log.write("%d subjects, %d rows in imputed growth data\n" %
                      (dbs.ID.unique().size, dbs.shape[0]))
            self.log.write("%d subjects, %d rows in cohort file data\n" %
                      (self.df.ID.unique().size, self.df.shape[0]))
            self.log.write("%d subjects, %d rows in merged data\n" %
                      (dx.ID.unique().size, dx.shape[0]))

        # Merge in imputed puberty variables
        #for vn in "Breast_Stage_Use_Z", "log2T_use_Z":
        #    dd = pd.read_csv(os.path.join("imputed_data_puberty", "%s_imp_%d.csv" % (vn, self.ix)))
        #    dd = pd.merge(dx, dd, left_on=("ID", "Age"), right_on=("ID", "Age"), how="left")
        #    ix = pd.notnull(dd[vn + "_x"])
        #    iy = pd.notnull(dd[vn + "_y"])
        #    dd[vn] = np.nan
        #    dd.loc[ix, vn] = dd.loc[ix, vn + "_x"]
        #    dd.loc[iy, vn] = dd.loc[iy, vn + "_y"]
        #    dx = dd.drop([vn + "_x", vn + "_y"], axis=1)
        #    if self.ix == 0:
        #        self.log.write("%d subjects, %d rows after merging %s\n" %
        #                  (dx.ID.unique().size, dx.shape[0], vn))

        # Code testosterone for females, breast stage for males as zero
        #dx.loc[dx.Female == 1, "log2T_use_Z"] = 0
        #dx.loc[dx.Female == 0, "Breast_Stage_Use_Z"] = 0

        dx = dx.loc[pd.notnull(dx[self.bp_var]), :]
        if self.ix == 0:
            self.log.write("%d blood pressure measures prior to final merge\n" % dx.shape[0])
            mv = pd.isnull(dx).sum(0)
            self.log.write("Missing values prior to final merge:\n%s\n" % mv.to_string())

        self.data = dx.dropna()

        if self.ix == 0:
            self.log.write("%d subjects, %d rows in final merged data\n" %
                      (self.data.ID.unique().size, self.data.shape[0]))
            dt = self.data[["ID", self.bp_var]].copy()
            nf = dt.groupby("ID")[self.bp_var].size()
            nf = pd.DataFrame(nf)
            nf = nf.reset_index()
            nf.columns = ["ID", "n_sbp"]
            nf.to_csv(os.path.join(self.bp_dir, "mixed", "stats", "%s_n.csv" % self.impvar),
                index=False)

        self.ix += 1
