import numpy as np
import pandas as pd
from scipy.stats import expon


def _get_distance_bw_contig(i, j, contigs):
    """
    Calculating distance between 2 contigs
    """
    distance_bw_contig = 0
    if i > j:
        i, j = j, i
    for k in range(i + 1, j):
        distance_bw_contig += len(contigs[k])

    return distance_bw_contig


def get_distance(pairs, contigs):
    """
    Calculating distances between all reads
    """
    L1 = np.array([len(contigs[int(pairs[i, 0])]) for i in range(len(pairs))])
    L2 = np.array([len(contigs[int(pairs[i, 2])]) for i in range(len(pairs))])

    S1 = pairs[:, 1] * (1 - pairs[:, 4]) + (L1 - pairs[:, 1] + 1) * (pairs[:, 4])
    S2 = (pairs[:, 3] + 1) * pairs[:, 5] + (L2 - pairs[:, 3]) * (1 - pairs[:, 5])

    S = S1 + S2 + pairs[:, 6]
    return S


def _distance_matrix(contigs):
    D = np.zeros((len(contigs), len(contigs)))

    for i in range(len(contigs)):
        for j in range(len(contigs)):
            D[i, j] = _get_distance_bw_contig(i, j, contigs)

    return D


def get_orientation(orientation, pairs, contigs):
    """
    Set up some orientation
    :param orientation: list of 1 for + or 0 for -
    """
    for i in range(len(contigs)):
        if contigs[i].o != orientation[i]:
            change_orientation(i, pairs, contigs)


def log_likelihood(pairs, contigs, P):
    """
    Calculating full log_lokelihood for our orientation
    """
    log_likelihood = P(get_distance(pairs, contigs)).sum()
    return log_likelihood


def change_orientation(number_changed_contig, pairs, contigs):
    """
    Changing orientation for one contig in pairs and in contigs
    """
    contigs[number_changed_contig].o = 1 if contigs[number_changed_contig].o == 0 else 0
    ind = contigs[number_changed_contig].reads_ind

    a = (pairs[ind][:, 0] == number_changed_contig)
    indx_1 = np.where(a)[0]
    indx_2 = np.where(~a)[0]

    pairs[ind[indx_1], 4] = -pairs[ind[indx_1], 4] + 1
    pairs[ind[indx_2], 5] = -pairs[ind[indx_2], 5] + 1


def change_log_likelihood(last_log_likelihood, number_changed_contig, pairs, contigs, P):
    """
    P(new_orientation) = P(old_orientation) + P(difference in orientation)
    """
    last_log_likelihood -= P(get_distance(pairs[contigs[number_changed_contig].reads_ind], contigs)).sum()
    change_orientation(number_changed_contig, pairs, contigs)
    new_lk = last_log_likelihood + P(get_distance(pairs[contigs[number_changed_contig].reads_ind], contigs)).sum()

    return new_lk


def simulation(max_len=1000, n_reads=10, n_contigs=2, p_distr=expon):
    """
    Creating data for simulation with our density
    :param max_len: max_len of genome
    :param n_reads: number of reads
    :param n_contigs: number of contigs
    :param p_distr: density
    """
    all_pos = []
    max_len = max_len - max_len % n_contigs
    len_cont = int(max_len / n_contigs)

    for i in range(n_reads):
        first = np.random.randint(0, max_len)
        second = first + int(p_distr.rvs() * 15)  # 4000
        if second < max_len:
            all_pos.append([first, second])

    # creating layout
    s = []
    for i in range(n_contigs):
        sym = "+"
        s.append("name{0}{1}".format(i, sym))
    with open("../data/simulation.layout.txt", "w") as f:
        f.writelines(",".join(s))

    # creating lens
    s = []
    for i in range(n_contigs):
        s.append("name{0}".format(i))
    z = [max_len // n_contigs for i in range(n_contigs)]
    df = pd.DataFrame([s, z]).T
    df.to_csv("../data/simulation.lens.tsv", sep="\t", index=False, header=False)

    # creating pairs
    d = ["\t".join(map(str, ["*",
                             s[all_pos[i][0] // len_cont],
                             all_pos[i][0] % len_cont,
                             s[all_pos[i][1] // len_cont],
                             all_pos[i][1] % len_cont,
                             "-", "-"])) + "\n"
         for i in range(len(all_pos))]
    with open("../data/simulation.pairs.txt", "w") as f:
        f.writelines(d)
