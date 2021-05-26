import numpy as np

from utils.tools import get_distance


def shuffle_ordering(pairs, contigs, P, n_iterations=10):
    """
    Random shuffling order if contigs
    """
    for _ in range(n_iterations):
        number_contig_1 = np.random.randint(0, len(contigs))
        number_contig_2 = np.random.randint(0, len(contigs))
        swap(number_contig_1, number_contig_2, pairs, contigs)


def swap(number_contig_1, number_contig_2, pairs, contigs):
    """
    Swapping 2 contigs and updating pairs
    """
    if number_contig_1 == number_contig_2:
        return
    if contigs[number_contig_1].pos > contigs[number_contig_2].pos:
        number_contig_1, number_contig_2 = number_contig_2, number_contig_1

    middle_length = 0
    for contig in contigs:
        if contigs[number_contig_1].pos < contig.pos < contigs[number_contig_2].pos:
            middle_length += contig.length

    for pair in pairs:
        if pair[0] == number_contig_1 and pair[2] == number_contig_2:
            pair[0], pair[2] = pair[2], pair[0]
            pair[1], pair[3] = pair[3], pair[1]
            pair[4], pair[5] = pair[5], pair[4]
        elif pair[0] == number_contig_1 and contigs[int(pair[2])].pos < contigs[number_contig_2].pos:
            pair[6] = middle_length - contigs[int(pair[2])].length - pair[6]
            pair[0], pair[2] = pair[2], pair[0]
            pair[1], pair[3] = pair[3], pair[1]
            pair[4], pair[5] = pair[5], pair[4]
        elif pair[0] == number_contig_1 and contigs[int(pair[2])].pos > contigs[number_contig_2].pos:
            pair[6] -= middle_length + contigs[number_contig_2].length
        elif pair[2] == number_contig_1:
            pair[6] += middle_length + contigs[number_contig_2].length
        elif pair[2] == number_contig_2 and contigs[int(pair[0])].pos > contigs[number_contig_1].pos:
            pair[6] = middle_length - contigs[int(pair[0])].length - pair[6]
            pair[0], pair[2] = pair[2], pair[0]
            pair[1], pair[3] = pair[3], pair[1]
            pair[4], pair[5] = pair[5], pair[4]
        elif pair[2] == number_contig_2 and contigs[int(pair[0])].pos < contigs[number_contig_1].pos:
            pair[6] -= middle_length + contigs[number_contig_1].length
        elif pair[0] == number_contig_2:
            pair[6] += middle_length + contigs[number_contig_1].length

    contigs[number_contig_1].pos, contigs[number_contig_2].pos = \
        contigs[number_contig_2].pos, contigs[number_contig_1].pos


def change_position(number_changed_contig, position_changed_contig, pairs, contigs):
    """
    Changing position of contig and updating pairs
    """
    if position_changed_contig == contigs[number_changed_contig].pos:
        return
    elif position_changed_contig < contigs[number_changed_contig].pos:
        while contigs[number_changed_contig].pos != position_changed_contig:
            for i in range(len(contigs)):
                if contigs[i].pos + 1 == contigs[number_changed_contig].pos:
                    swap(i, number_changed_contig, pairs, contigs)
                    break
    else:
        while contigs[number_changed_contig].pos != position_changed_contig:
            for i in range(len(contigs)):
                if contigs[i].pos - 1 == contigs[number_changed_contig].pos:
                    swap(i, number_changed_contig, pairs, contigs)
                    break


def swap_log_likelihood(last_log_likelihood, number_contig_1, number_contig_2, pairs, contigs, P):
    """
    P(new_orientation) = P(old_orientation) + P(difference in orientation)
    """
    reads_ind = list(set(contigs[number_contig_1].reads_ind).union(set(contigs[number_contig_2].reads_ind)))

    last_log_likelihood -= P(get_distance(pairs[reads_ind], contigs)).sum()
    swap(number_contig_1, number_contig_2, pairs, contigs)
    new_lk = last_log_likelihood
    new_lk += P(get_distance(pairs[reads_ind], contigs)).sum()

    return new_lk


def change_position_log_likelihood(last_log_likelihood, number_changed_contig, position_changed_contig, pairs, contigs,
                                   P):
    """
    P(new_orientation) = P(old_orientation) + P(difference in orientation)
    """
    reads_ind = set(contigs[number_changed_contig].reads_ind)

    left_pos = min(position_changed_contig, contigs[number_changed_contig].pos)
    right_pos = max(position_changed_contig, contigs[number_changed_contig].pos)

    for i in range(len(contigs)):
        if left_pos <= contigs[i].pos <= right_pos:
            reads_ind.union(set(contigs[i].reads_ind))
    reads_ind = list(reads_ind)

    last_log_likelihood -= P(get_distance(pairs[reads_ind], contigs)).sum()
    change_position(number_changed_contig, position_changed_contig, pairs, contigs)
    new_lk = last_log_likelihood
    new_lk += P(get_distance(pairs[reads_ind], contigs)).sum()

    return new_lk
