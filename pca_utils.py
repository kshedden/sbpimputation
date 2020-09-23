import numpy as np
import pandas as pd
import os

def do_pca(impvar, ndim, bp_dir, plot=False):

    # The variable names for the trajectory of childhood exposures
    vx = [(impvar + "%d") % a for a in range(1, 11)]

    # Second derivatives for smoothing the PC's
    p = len(vx)
    fx = np.zeros((p-2, p))
    for i in range(p-2):
        fx[i, i:i+3] = [1, -2, 1]

    # Get the moments of the imputed data
    mn, cov = 0, 0
    for k in range(20):
        di = pd.read_csv(os.path.join("imputed_data", "%s_imp_%d.csv" % (impvar, k)))
        di["ID"] = di["ID"].astype(np.int)
        dd = di[vx]
        mn += dd.mean(0)
        cov += np.cov(dd.T)
    mn /= 20
    cov /= 20

    # Functional PCA
    xsd = np.sqrt(np.diag(cov))
    cor = cov / np.outer(xsd, xsd)
    cm = cor - 0.2*np.dot(fx.T, fx)
    b, proj = np.linalg.eig(cm)
    ii = np.argsort(b)[::-1]
    b = b[ii]
    proj = proj[:, ii]
    b = b[0:ndim]
    proj = proj[:, 0:ndim]
    proj = proj / xsd[:, None]
    if np.any(b < 0):
        1/0

    # Plot basis functions
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    d = "mixed"
    pdf = PdfPages(os.path.join(bp_dir, d, "dim_%d" % ndim, "%s_basis.pdf" % impvar))
    plt.clf()
    plt.title("%s basis functions" % impvar)
    plt.grid(True)
    for k in range(ndim):
        plt.plot(np.arange(1, proj.shape[0]+1), proj[:, k], '-', color='orange', lw=4)
    plt.xlabel("Age", size=15)
    pdf.savefig()
    pdf.close()

    # The names of the PC score variables
    vb = ["%s_pc%d" % (impvar, j) for j in range(proj.shape[1])]

    return mn, cov, proj, vb, vx
