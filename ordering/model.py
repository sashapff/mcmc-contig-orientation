import numpy as np
import random
from tqdm import tqdm
from ordering.tools import get_ordering, change_position, change_position_log_likelihood
from tools.tools import log_likelihood


def MCMC(pairs, contigs, P, number_it=500):
    """
    Changing orientation for better likelihood
    :param pairs: array of reads
    :param contigs: list of contigs
    :param P: density
    :param number_it: number of iterations
    :param n_chains: number of changing for one iterations
    :return None, change on place
    """
    lk_old = log_likelihood(pairs, contigs, P)
    new_contigs = [contig.pos for contig in contigs]

    log_likelihood_arr = []

    for _ in tqdm(range(number_it)):
        actual_shifts = 0
        while actual_shifts == 0:
            number_changed_contig = np.random.randint(0, len(contigs))
            shifts_changed_contig = np.random.randint(-5, 5)
            if shifts_changed_contig >= 0:
                shifts_changed_contig += 1

            new_position = min(len(contigs) - 1, max(0, contigs[number_changed_contig].pos + shifts_changed_contig))
            actual_shifts = contigs[number_changed_contig].pos - new_position

        lk_new = change_position_log_likelihood(lk_old, number_changed_contig, new_position, pairs, contigs, P)
        if actual_shifts == 0:
            assert np.abs(lk_old - lk_new) < 1

        if random.random() > np.exp(lk_new - lk_old):
            # Decline
            assert 0 <= contigs[number_changed_contig].pos + actual_shifts < len(contigs)
            change_position(number_changed_contig, contigs[number_changed_contig].pos + actual_shifts, pairs, contigs)
        else:
            # Accept
            lk_old = lk_new
            new_contigs = [contig.pos for contig in contigs]

        log_likelihood_arr.append(lk_old)
        print([contig.pos for contig in contigs], lk_old)

    # assert new_contigs == [contig.pos for contig in contigs]
    # get_ordering(new_contigs, pairs, contigs)

    return log_likelihood_arr
