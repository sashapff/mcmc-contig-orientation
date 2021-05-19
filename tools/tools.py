import numpy as np
from scipy.stats import expon, gamma
from tqdm import tqdm
import pandas as pd


def get_longest_contig(longest_contig_pairs_arr, longest_contig_arr):
    j = 0
    longest_contig_length = 0
    for i in range(len(longest_contig_arr)):
        if longest_contig_arr[i].length > longest_contig_length:
            longest_contig_length = longest_contig_arr[i].length
            j = i

    return longest_contig_pairs_arr[j]


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


def get_distance_one_contig(filtered_pairs, left_len):
    """
    Calculating distances between all reads
    """
    L1 = left_len

    S1 = (L1 - filtered_pairs[:, 1] + 1)
    S2 = filtered_pairs[:, 3] + 1

    S = S1 + S2 + filtered_pairs[:, 6]
    return S


def _distance_matrix(contigs):
    D = np.zeros((len(contigs), len(contigs)))

    for i in range(len(contigs)):
        for j in range(len(contigs)):
            D[i, j] = _get_distance_bw_contig(i, j, contigs)

    return D


def log_likelihood(pairs, contigs, P):
    """
    Calculating full log_likelihood for our orientation
    """
    return P(get_distance(pairs, contigs)).sum()


def filter_pairs(pairs, contig_id, left, right):
    indx = (pairs[:, 0] == contig_id) & (pairs[:, 2] == contig_id) & (pairs[:, 1] < left) & (pairs[:, 3] >= right)
    # indx = (pairs[:, 0] == contig_id) & (pairs[:, 2] == contig_id) & (pairs[:, 3] - pairs[:, 1] >= d)
    return pairs[indx]


def simulation(max_len=1000, n_reads=10, n_contigs=2, p_distr=expon, output_path='../data'):
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

    for i in tqdm(range(n_reads)):
        first = np.random.randint(0, max_len)
        second = first + int(p_distr.rvs() * 150)  # 4000
        if second < max_len:
            all_pos.append([first, second])

    # creating layout
    s = []
    for i in range(n_contigs):
        sym = "+"
        s.append("name{0}{1}".format(i, sym))
    with open(output_path + "simulation.layout.txt", "w") as f:
        f.write('* ')
        f.writelines(",".join(s))

    # creating lens
    s = []
    for i in range(n_contigs):
        s.append("name{0}".format(i))
    z = [max_len // n_contigs for i in range(n_contigs)]
    df = pd.DataFrame([s, z]).T
    df.to_csv(output_path + "simulation.lens.tsv", sep="\t", index=False, header=False)

    # creating pairs
    d = ["\t".join(map(str, ["*",
                             s[all_pos[i][0] // len_cont],
                             all_pos[i][0] % len_cont,
                             s[all_pos[i][1] // len_cont],
                             all_pos[i][1] % len_cont,
                             "-", "-"])) + "\n"
         for i in range(len(all_pos))]
    with open(output_path + "simulation.pairs.txt", "w") as f:
        f.writelines(d)
