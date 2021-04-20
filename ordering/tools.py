import numpy as np

from tools.tools import get_distance, _get_distance_bw_contig, log_likelihood


def shuffle_ordering(pairs, contigs, P, n_iterations=10):
    """
    TODO
    """
    for _ in range(n_iterations):
        number_contig_1 = np.random.randint(0, len(contigs))
        number_contig_2 = np.random.randint(0, len(contigs))
        swap(number_contig_1, number_contig_2, pairs, contigs)
        # swap_neigh(number_contig_1, pairs, contigs)

    # for pos in range(len(contigs)):
    #     for i in range(len(contigs)):
    #         if ordering[i] == pos and contigs[i].pos != ordering[i]:
    #             change_position(i, pos, pairs, contigs)
    #             break


def swap(number_contig_1, number_contig_2, pairs, contigs):
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


def swap_neigh(number_contig_1, pairs, contigs):
    number_contig_2 = number_contig_1

    if contigs[number_contig_1].pos != len(contigs) - 1:
        for (i, contig) in enumerate(contigs):
            if contig.pos == contigs[number_contig_1].pos + 1:
                number_contig_2 = i
    else:
        for (i, contig) in enumerate(contigs):
            if contig.pos == contigs[number_contig_1].pos - 1:
                number_contig_2 = i
        number_contig_1, number_contig_2 = number_contig_2, number_contig_1

    for pair in pairs:
        if pair[0] == number_contig_1 and pair[2] == number_contig_2:
            pair[0], pair[2] = pair[2], pair[0]
            pair[1], pair[3] = pair[3], pair[1]
            pair[4], pair[5] = pair[5], pair[4]
        # elif pair[0] == number_contig_1 and contigs[int(pair[2])].pos < contigs[number_contig_2].pos:
        #     pair[6] = middle_length - contigs[int(pair[2])].length - pair[6]
        #     pair[0], pair[2] = pair[2], pair[0]
        #     pair[1], pair[3] = pair[3], pair[1]
        #     pair[4], pair[5] = pair[5], pair[4]
        elif pair[0] == number_contig_1 and contigs[int(pair[2])].pos > contigs[number_contig_2].pos:
            pair[6] -= contigs[number_contig_2].length
        elif pair[2] == number_contig_1:
            pair[6] += contigs[number_contig_2].length
        # elif pair[2] == number_contig_2 and contigs[int(pair[0])].pos > contigs[number_contig_1].pos:
        #     pair[6] = middle_length - contigs[int(pair[0])].length - pair[6]
        #     pair[0], pair[2] = pair[2], pair[0]
        #     pair[1], pair[3] = pair[3], pair[1]
        #     pair[4], pair[5] = pair[5], pair[4]
        elif pair[2] == number_contig_2 and contigs[int(pair[0])].pos < contigs[number_contig_1].pos:
            pair[6] -= contigs[number_contig_1].length
        elif pair[0] == number_contig_2:
            pair[6] += contigs[number_contig_1].length

    contigs[number_contig_1].pos, contigs[number_contig_2].pos = \
        contigs[number_contig_2].pos, contigs[number_contig_1].pos


def change_position(number_changed_contig, position_changed_contig, pairs, contigs):
    """
    TODO
    """
    if position_changed_contig < contigs[number_changed_contig].pos:
        previous_pos = contigs[number_changed_contig].pos
        shift_length = 0
        for i in range(len(contigs)):
            if position_changed_contig <= contigs[i].pos < previous_pos:
                contigs[i].pos += 1

                ind = contigs[i].reads_ind
                a = (pairs[ind][:, 0] == i)
                indx_1 = np.where(a)[0]
                indx_2 = np.where(~a)[0]

                pairs[ind[indx_1], 6] -= contigs[number_changed_contig].length
                pairs[ind[indx_2], 6] += contigs[number_changed_contig].length

                # delta -= contigs[number_changed_contig].length * len(pairs[ind[indx_1], 6])
                # delta += contigs[number_changed_contig].length * len(pairs[ind[indx_2], 6])

                shift_length += contigs[i].length

        ind = contigs[number_changed_contig].reads_ind
        a = (pairs[ind][:, 0] == number_changed_contig)
        indx_1 = np.where(a)[0]
        indx_2 = np.where(~a)[0]

        pairs[ind[indx_1], 6] += shift_length
        pairs[ind[indx_2], 6] -= shift_length

        # delta += shift_length * len(pairs[ind[indx_1], 6])
        # delta -= shift_length * len(pairs[ind[indx_2], 6])

        for i in range(len(contigs)):
            if position_changed_contig < contigs[i].pos <= previous_pos and i != number_changed_contig:
                ind = contigs[i].reads_ind
                a = (pairs[ind][:, 0] == i)
                indx0 = np.where(a)[0]
                b = (pairs[ind[indx0]][:, 2] == number_changed_contig)
                indx = np.where(b)[0]

                # delta -= 2 * len(pairs[ind[indx], 6])

                pairs[ind[indx], 6] = -pairs[ind[indx], 6]

                pairs[ind[indx], 0], pairs[ind[indx], 2] = pairs[ind[indx], 2], pairs[ind[indx], 0]
                pairs[ind[indx], 1], pairs[ind[indx], 3] = pairs[ind[indx], 3], pairs[ind[indx], 1]
                pairs[ind[indx], 4], pairs[ind[indx], 5] = pairs[ind[indx], 5], pairs[ind[indx], 4]

        contigs[number_changed_contig].pos = position_changed_contig
    else:
        previous_pos = contigs[number_changed_contig].pos
        shift_length = 0
        for i in range(len(contigs)):
            if position_changed_contig >= contigs[i].pos > previous_pos:
                contigs[i].pos -= 1

                ind = contigs[i].reads_ind
                a = (pairs[ind][:, 0] == i)
                indx_1 = np.where(a)[0]
                indx_2 = np.where(~a)[0]

                pairs[ind[indx_1], 6] += contigs[number_changed_contig].length
                pairs[ind[indx_2], 6] -= contigs[number_changed_contig].length

                # delta += contigs[number_changed_contig].length * len(pairs[ind[indx_1], 6])
                # delta -= contigs[number_changed_contig].length * len(pairs[ind[indx_2], 6])

                shift_length += contigs[i].length

        ind = contigs[number_changed_contig].reads_ind
        a = (pairs[ind][:, 0] == number_changed_contig)
        indx_1 = np.where(a)[0]
        indx_2 = np.where(~a)[0]

        pairs[ind[indx_1], 6] -= shift_length
        pairs[ind[indx_2], 6] += shift_length

        # delta -= shift_length * len(pairs[ind[indx_1], 6])
        # delta += shift_length * len(pairs[ind[indx_2], 6])

        for i in range(len(contigs)):
            if position_changed_contig > contigs[i].pos >= previous_pos and i != number_changed_contig:
                ind = contigs[i].reads_ind
                a = (pairs[ind][:, 2] == i)
                indx0 = np.where(a)[0]
                b = (pairs[ind[indx0]][:, 0] == number_changed_contig)
                indx = np.where(b)[0]

                # delta -= 2 * len(pairs[ind[indx], 6])

                pairs[ind[indx], 6] = -pairs[ind[indx], 6]

                pairs[ind[indx], 0], pairs[ind[indx], 2] = pairs[ind[indx], 2], pairs[ind[indx], 0]
                pairs[ind[indx], 1], pairs[ind[indx], 3] = pairs[ind[indx], 3], pairs[ind[indx], 1]
                pairs[ind[indx], 4], pairs[ind[indx], 5] = pairs[ind[indx], 5], pairs[ind[indx], 4]

        contigs[number_changed_contig].pos = position_changed_contig

    # return delta


def swap_log_likelihood(last_log_likelihood, number_contig_1, number_contig_2, pairs, contigs, P):
    """
    P(new_orientation) = P(old_orientation) + P(difference in orientation)
    """
    last_log_likelihood -= P(get_distance(pairs[contigs[number_contig_1].reads_ind], contigs)).sum()
    last_log_likelihood -= P(get_distance(pairs[contigs[number_contig_2].reads_ind], contigs)).sum()

    swap(number_contig_1, number_contig_2, pairs, contigs)

    new_lk = last_log_likelihood
    new_lk += P(get_distance(pairs[contigs[number_contig_1].reads_ind], contigs)).sum()
    new_lk += P(get_distance(pairs[contigs[number_contig_2].reads_ind], contigs)).sum()

    return new_lk


def change_position_log_likelihood(last_log_likelihood, number_changed_contig, position_changed_contig, pairs, contigs,
                                   P):
    """
    P(new_orientation) = P(old_orientation) + P(difference in orientation)
    """
    # if contigs[number_changed_contig].pos < position_changed_contig:
    #     left = contigs[number_changed_contig].pos
    #     right = position_changed_contig + 1
    # else:
    #     left = position_changed_contig
    #     right = contigs[number_changed_contig].pos + 1
    # changed_contigs = set()
    # for i in range(left, right):
    #     changed_contigs = changed_contigs.union(set(contigs[i].reads_ind))
    # changed_contigs = np.array(list(changed_contigs))
    # last_log_likelihood -= P(get_distance(pairs[changed_contigs], contigs)).sum()
    # old_pos = contigs[number_changed_contig].pos
    # delta = change_position(number_changed_contig, position_changed_contig, pairs, contigs, 0)
    # print(delta)
    # delta = change_position(number_changed_contig, old_pos, pairs, contigs, 0)
    # print(delta)
    # new_lk = last_log_likelihood + P(get_distance(pairs[changed_contigs], contigs)).sum()
    # return new_lk
