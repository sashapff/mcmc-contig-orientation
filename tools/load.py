import numpy as np
import pandas as pd
from tools.tools import _distance_matrix
from tqdm import tqdm


class Contig():
    def __init__(self, name, length, orientation, position):
        """
        class of contig
        :param name: name
        :param length: length of contig
        :param orientation: orientation of contig
        TODO
        """
        self.name = name
        self.length = length
        self.o = orientation
        self.pos = position
        self.reads_ind = None

    def __len__(self):
        return self.length

    def define_reads_in_contig(self, pairs, index):
        self.reads_ind = np.argwhere(((pairs[:, 0] == index) | (pairs[:, 2] == index))
                                     & (pairs[:, 0] != pairs[:, 2])).T[0]


def _clear_layout(path_layout, path_lens, min_len=10e+5):
    """
    Reading layout and deleting contigs which less than min_len
    """
    with open(path_layout, "r") as f:
        s_lines = f.read().splitlines()

    length = pd.read_csv(path_lens, sep="\t", header=None)
    print("COULD ^C")

    s = s_lines[0]
    for line in s_lines:
        if len(line) > len(s):
            s = line

    s = (s.split(" ")[1]).split(",")

    i = 0
    while i != len(s):
        if int(length.loc[length[0] == s[i][:-1], 1]) < min_len:
            s.remove(s[i])
        else:
            i += 1
    return s

def check_reads(path_pairs):
    with open(path_pairs, "r") as f:
        s_lines = f.read().splitlines()

    cnt = 0
    for line in tqdm(s_lines):
        line = line.split('\t')
        if line[1] != line[3]:
            cnt += 1
    return cnt


def get_contigs_and_pairs(path_layout, path_lens, path_pairs, long_contig=False, all_contigs=False, min_len=100_000):
    """
    Reading all data from file in essential format
    :param path_layout: the path to the start layout
    :param path_lens: the path to the length of contigs
    :param path_pairs: all reads
    :param long_contig: bool, return reads from longest contig separately?
    :param min_len: min length of contigs we will deal with
    :param all_contigs: return all reads in same contigs?
    :return: required reads, list of contigs, dict{"name": index contig in list},
             long_contig if True, all_contigs if True
    """
    length = pd.read_csv(path_lens, sep="\t", header=None)
    answer = []

    s = _clear_layout(path_layout, path_lens, min_len=min_len)

    print("Reading pairs...")
    pairs = pd.read_csv(path_pairs, "\t", names=["name", "X1", "P1", "X2", "P2", "orientation1", "orientation2"])

    # list of contigs with right layout
    contigs = [Contig(s[i][:-1], int(length[length[0] == s[i][:-1]][1]), int(s[i][-1] == "+"), i) for i in
               range(0, len(s))]

    # dict "name contig" -> id contig in contigs
    id_contig = {contig.name: i for i, contig in enumerate(contigs)}

    print("Cleaning pairs...")
    # to remove pairs inside of contigs and which in contigs are not suitable for the length
    if long_contig:
        longest_contig = max(contigs, key=len).name
        indx = (pairs["X1"] == longest_contig) & (pairs["X2"] == longest_contig)
        longest_pairs_numpy = np.zeros((indx.sum(), 2))

        longest_pairs_numpy[:, 0] = pairs[indx]["P1"].to_numpy()
        longest_pairs_numpy[:, 1] = pairs[indx]["P2"].to_numpy()

        del indx
        answer.append(longest_pairs_numpy)
        answer.append(max(contigs, key=len))

    if all_contigs:
        contigs_list = []
        for contig in contigs:
            indx = (pairs["X1"] == contig.name) & (pairs["X2"] == contig.name)
            pairs_numpy = np.zeros((indx.sum(), 2))

            pairs_numpy[:, 0] = pairs[indx]["P1"].to_numpy()
            pairs_numpy[:, 1] = pairs[indx]["P2"].to_numpy()
            del indx
            contigs_list.append(pairs_numpy)
        answer.append(contigs_list)

    pairs = pairs[pairs["X1"] != pairs["X2"]]
    pairs = pairs[
        pairs["X1"].apply(lambda x: x in id_contig) & pairs["X2"].apply(lambda x: x in id_contig)].reset_index(
        drop=True)

    print("Counting distance...")
    D = _distance_matrix(contigs)

    np_pairs = np.zeros((len(pairs), 7))
    np_pairs[:, 0] = pairs["X1"].apply(lambda x: id_contig[x]).to_numpy()
    np_pairs[:, 1] = pairs["P1"].to_numpy()
    np_pairs[:, 2] = pairs["X2"].apply(lambda x: id_contig[x]).to_numpy()
    np_pairs[:, 3] = pairs["P2"].to_numpy()
    np_pairs[:, 4] = pairs["X1"].apply(lambda x: contigs[id_contig[x]].o).to_numpy()
    np_pairs[:, 5] = pairs["X2"].apply(lambda x: contigs[id_contig[x]].o).to_numpy()
    np_pairs[:, 6] = np.array([D[int(np_pairs[i, 0]), int(np_pairs[i, 2])] for i in range(len(pairs))])

    # id_contig1 < id_contig2
    ind = (np_pairs[:, 0] >= np_pairs[:, 2])
    np_pairs[ind, 0], np_pairs[ind, 2] = np_pairs[ind, 2], np_pairs[ind, 0]
    np_pairs[ind, 1], np_pairs[ind, 3] = np_pairs[ind, 3], np_pairs[ind, 1]
    np_pairs[ind, 4], np_pairs[ind, 5] = np_pairs[ind, 5], np_pairs[ind, 4]

    print("Reading reads for contigs...")
    for index, contig in enumerate(contigs):
        contig.define_reads_in_contig(np_pairs, index)

    del pairs
    print("Done!")

    return [np_pairs, contigs, id_contig] + answer
