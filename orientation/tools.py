import numpy as np
from tools.tools import get_distance


def get_orientation(orientation, pairs, contigs):
    """
    Set up some orientation
    :param orientation: list of 1 for + or 0 for -
    """
    for i in range(len(contigs)):
        if contigs[i].o != orientation[i]:
            change_orientation(i, pairs, contigs)


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


def change_orientation_log_likelihood(last_log_likelihood, number_changed_contig, pairs, contigs, P):
    """
    P(new_orientation) = P(old_orientation) + P(difference in orientation)
    """
    last_log_likelihood -= P(get_distance(pairs[contigs[number_changed_contig].reads_ind], contigs)).sum()
    change_orientation(number_changed_contig, pairs, contigs)
    new_lk = last_log_likelihood + P(get_distance(pairs[contigs[number_changed_contig].reads_ind], contigs)).sum()

    return new_lk
