import numpy as np

from tools.tools import get_distance, _get_distance_bw_contig


def get_ordering(ordering, pairs, contigs):
    """
    TODO
    """
    for i in range(len(contigs)):
        for j in range(i + 1, len(contigs)):
            if ordering[i] > ordering[j]:
                swap(i, j, pairs, contigs)

    # for pos in range(len(contigs)):
    #     for i in range(len(contigs)):
    #         if ordering[i] == pos and contigs[i].pos != ordering[i]:
    #             change_position(i, pos, pairs, contigs)
    #             break


def swap(number_contig_1, number_contig_2, pairs, contigs):
    if contigs[number_contig_1].pos > contigs[number_contig_2].pos:
        number_contig_1, number_contig_2 = number_contig_2, number_contig_1

    # middle_contigs = []
    middle_length = _get_distance_bw_contig(number_contig_1, number_contig_2, contigs)
    # for (i, contig) in enumerate(contigs):
    #     if contigs[number_contig_1].pos < contig.pos < contigs[number_contig_2].pos:
    #         # middle_contigs.append(i)
    #         middle_length += contig.length
    # middle_contigs = np.array(middle_contigs)

    # ind_1 = contigs[number_contig_1].reads_ind
    # ind_2 = contigs[number_contig_2].reads_ind
    #
    # indx_1_left = np.where((pairs[ind_1][:, 2] == number_contig_1))[0]
    # indx_1_right = np.where((pairs[ind_1][:, 0] == number_contig_1))[0]
    # indx_2_left = np.where((pairs[ind_2][:, 2] == number_contig_2))[0]
    # indx_2_right = np.where((pairs[ind_2][:, 0] == number_contig_2))[0]
    # #
    # # indx_1_2 = np.where((pairs[ind_1][:, 0] == number_contig_1) and (pairs[ind_1][:, 2] == number_contig_2))[0]
    # #
    # pairs[ind_1[indx_1_left], 6] += middle_length
    # pairs[ind_2[indx_2_right], 6] += middle_length
    #
    # pairs[ind_1[indx_1_right], 6] -= middle_length
    # pairs[ind_2[indx_2_left], 6] -= middle_length

    for pair in pairs:
        if pair[0] == number_contig_1 and pair[2] == number_contig_2:
            pass
        if pair[0] == number_contig_1 and contigs[pair[2]].pos < contigs[number_contig_2].pos:
            pair[6] = -(pair[6] + contigs[pair[2]].length)
        if pair[0] == number_contig_1 and contigs[pair[2]].pos > contigs[number_contig_2].pos:
            pair[6] -= middle_length
        if pair[2] == number_contig_1:
            pair[6] += middle_length
        if pair[2] == number_contig_2 and contigs[pair[0]].pos > contigs[number_contig_1].pos:
            pair[6] = -(pair[6] + contigs[pair[0]].length)
        if pair[2] == number_contig_2 and contigs[pair[0]].pos < contigs[number_contig_1].pos:
            pair[6] -= middle_length
        if pair[0] == number_contig_2:
            pair[6] += middle_length


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
