import numpy as np
from sklearn.neighbors.kde import KernelDensity
from scipy.optimize import curve_fit


def density(reads, K0=3000, K1=100_000, K2=400_000, kde_method="linear"):
    """
    Estimating density of distances between peace of reads
    :param reads: reads[i,0] - pisition first peace of read; reads[i,1] - pisition second peace of read
    :param K0: end of poly approximation
    :param K1: end of first log approximation
    :param K2: end of second log approximation
    :param kde_method: type of kernel in KernelDensity
    :return:
    """
    distances = np.abs(reads[:, 0] - reads[:, 1])

    kde = KernelDensity(kernel=kde_method, bandwidth=200).fit(distances.reshape(-1, 1))
    f = lambda x: kde.score_samples(x.reshape(-1, 1))

    x1 = np.logspace(np.log10(K0 - 1000) , np.log10(K1), 500)
    p = lambda x, a, b: a + b*np.log(x)
    param1, cov = curve_fit(p, x1, f(x1))

    # proximal
    degree = 30
    x0 = np.logspace(0, np.log10(K0 + 1000), 500)
    param0 = np.polyfit(x0, f(x0), degree)

    x2 = np.logspace(np.log10(K1 - 1000) , np.log10(K2), 500)
    p = lambda x, a, b: a + b*np.log(x)
    param2, cov = curve_fit(p, x2, f(x2))

    P = (lambda x: np.where(x < K0, np.poly1d(param0)(x),
                   np.where(x < K1, param1[0] + param1[1] * np.log(x),
                   np.where(x < K2, param2[0] + param2[1]*np.log(x),
                                    param2[0] + param2[1]*np.log(K2)))))

    return P, f


def toy_density(longest_contig):
    """
    Function for estimating simulation density
    """
    distances = np.abs(longest_contig[:, 0] - longest_contig[:, 1])
    Lambda = 1 / distances.mean()
    P = lambda x: Lambda * np.exp(-Lambda * x)

    return P, P
