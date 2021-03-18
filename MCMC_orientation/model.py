import numpy as np
import random
from tqdm import tqdm
from tools import log_likelihood, change_log_likelihood, change_orientation, get_orientation


def MCMC(pairs, contigs, P, number_it=500, n_chains=1):
    """
    Changing orientation for better likelihood
    :param pairs: array of reads
    :param contigs: list of contigs
    :param P: density
    :param number_it: number of itterations
    :param n_chains: number of changing for one itterations
    :return None, change on place
    """
    lk_old = log_likelihood(pairs, contigs, P)
    new_contigs = [contig.o for contig in contigs]

    for _ in tqdm(range(number_it)):
        numbers_changed_contig = np.random.randint(0, len(contigs), n_chains)

        for number in numbers_changed_contig:
            lk_new = change_log_likelihood(lk_old, number, pairs, contigs, P)

        if random.random() > np.exp(lk_new - lk_old):
            for number in numbers_changed_contig:
                change_orientation(number, pairs, contigs)
        else:
            lk_old = lk_new
            new_contigs = [contig.o for contig in contigs]

    get_orientation(new_contigs, pairs, contigs)

    return new_contigs
