import random

import numpy as np
from tqdm import tqdm

from ordering.tools import change_position, swap, swap_log_likelihood
from utils.tools import log_likelihood


def MCMC(pairs, contigs, P, number_it=500, strategy='swap'):
    """
    Changing orientation for better likelihood
    :param pairs: array of reads
    :param contigs: list of contigs
    :param P: density
    :param number_it: number of iterations
    :param strategy: strategy of random order changes
    :return None, change on place
    """
    lk_old = log_likelihood(pairs, contigs, P)

    print('Shuffled contigs:', lk_old, [contig.pos for contig in contigs])

    log_likelihood_arr = []

    for _ in tqdm(range(number_it)):
        number_contig_1 = np.random.randint(0, len(contigs))
        number_contig_2 = np.random.randint(0, len(contigs))

        last_contig_pos = contigs[number_contig_1].pos

        if strategy == 'swap':
            lk_new = swap_log_likelihood(lk_old, number_contig_1, number_contig_2, pairs, contigs, P)
        else:
            change_position(number_contig_1, number_contig_2, pairs, contigs)
            lk_new = log_likelihood(pairs, contigs, P)

        if random.random() > np.exp(lk_new - lk_old):
            # Decline
            if strategy == 'swap':
                swap(number_contig_1, number_contig_2, pairs, contigs)
            else:
                change_position(number_contig_1, last_contig_pos, pairs, contigs)
        else:
            # Accept
            lk_old = lk_new

        log_likelihood_arr.append(lk_old)

    return log_likelihood_arr
