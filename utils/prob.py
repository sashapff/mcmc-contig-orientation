import numpy as np
from scipy.integrate import quad
from scipy.optimize import curve_fit
from sklearn.neighbors.kde import KernelDensity


def normalize(P, a=0, b=np.inf):
    """
    Normalize the tail of the distribution to the area value under plot 1.
    """
    quad_value, _ = quad(lambda x: np.exp(P(x)), a, b)
    coeff = np.log(quad_value)
    P_norm = lambda x: P(x) - coeff
    return P_norm


def density(reads, K0=3000, K1=100_000, kde_method="linear"):
    """
    Estimating density of distances between peace of reads
    :param reads: reads[i,0] - position first peace of read; reads[i,1] - position second peace of read
    :param K0: end of poly approximation
    :param K1: end of first log approximation
    :param K2: end of second log approximation
    :param kde_method: type of kernel in KernelDensity
    :return:
    """
    distances = np.abs(reads[:, 0] - reads[:, 1])

    kde = KernelDensity(kernel=kde_method, bandwidth=200).fit(distances.reshape(-1, 1))
    f = lambda x: kde.score_samples(x.reshape(-1, 1))

    # proximal
    degree = 30
    x0 = np.logspace(0, np.log10(K0 + 1000), 500)
    param0 = np.polyfit(x0, f(x0), degree)

    x1 = np.logspace(np.log10(K0 - 1000), np.log10(K1), 500)
    p = lambda x, a, b: a + b * np.log(x)
    param1, cov = curve_fit(p, x1, f(x1))

    P = (lambda x: np.where(x < K0, np.poly1d(param0)(x),
                            np.where(x < K1, param1[0] + param1[1] * np.log(x),
                                     param1[0] + param1[1] * np.log(x))))

    return P, f


def toy_density(reads):
    """
    Function for estimating simulation density
    """
    distances = np.abs(reads[:, 0] - reads[:, 1])
    Lambda = 1 / distances.mean()
    P = lambda x: np.log(Lambda) - Lambda * x

    return P, P


def simulate_density(pairs):
    """
    Function for estimating simulation density
    """
    distances = np.abs(pairs[:, 1] - pairs[:, 3])
    Lambda = 1 / distances.mean()
    P = lambda x: np.log(Lambda) - Lambda * x
    return P, P


def estimate_density(reads, K0=3000, K1=100_000, kde_method="linear"):
    """
    Estimating density of distances between peace of reads
    :param reads: reads[i,0] - position first peace of read; reads[i,1] - position second peace of read
    :param K0: end of poly approximation
    :param K1: end of first log approximation
    :param kde_method: type of kernel in KernelDensity
    :return:
    """
    distances = np.abs(reads[:, 1] - reads[:, 3])

    kde = KernelDensity(kernel=kde_method, bandwidth=200).fit(distances.reshape(-1, 1))
    f = lambda x: kde.score_samples(x.reshape(-1, 1))

    # proximal
    degree = 30
    x0 = np.logspace(0, np.log10(K0 + 1000), 500)
    param0 = np.polyfit(x0, f(x0), degree)

    x1 = np.logspace(np.log10(K0 - 1000), np.log10(K1), 500)
    p = lambda x, a, b: a + b * np.log(x)
    param1, cov = curve_fit(p, x1, f(x1))

    P = (lambda x: np.where(x < K0, np.poly1d(param0)(x),
                            np.where(x < K1, param1[0] + param1[1] * np.log(x),
                                     param1[0] + param1[1] * np.log(x))))

    return P, f


def destiny_b(longest_contig_b, bins, contigs, resolution="100_000"):
    value, first_bins, second_bins = longest_contig_b

    long_id = np.argmax([len(contig) for contig in contigs])
    seq = bins.seq[bins.con_in_seq[long_id][0]: bins.con_in_seq[long_id][1]]
    long_len = len(contigs[long_id]) // int(resolution)

    means = np.zeros(len(seq))
    for i in range(len(seq)):
        means[i] = value[(second_bins - first_bins) == i].sum() / (long_len - i)

    return means[1:]
