import numpy as np
import random
from tqdm import tqdm
from ordering.tools import get_ordering, change_position, change_position_log_likelihood
from utils.tools import log_likelihood


def MCMC(pairs, contigs, P, number_it=500, n_chains=1):
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

    for _ in tqdm(range(number_it)):
        numbers_changed_contig = np.random.randint(0, len(contigs), n_chains)
        shifts_changed_contig = np.random.randint(-5, 5, n_chains)

        for (i, number) in enumerate(numbers_changed_contig):
            shift = shifts_changed_contig[i]
            lk_new = change_position_log_likelihood(lk_old, number,
                                                    min(len(contigs) - 1, max(0, contigs[number].pos + shift)), pairs,
                                                    contigs, P)

        if random.random() > np.exp(lk_new - lk_old):
            for (i, number) in enumerate(numbers_changed_contig):
                shift = shifts_changed_contig[i]
                change_position(number, min(len(contigs) - 1, max(0, contigs[number].pos + shift)), pairs, contigs)
        else:
            lk_old = lk_new
            new_contigs = [contig.pos for contig in contigs]

    get_ordering(new_contigs, pairs, contigs)

    return new_contigs
