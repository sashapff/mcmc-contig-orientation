import numpy as np


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
    return pairs[indx]
