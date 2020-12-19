from .load import get_contigs_and_pairs
from .prob import density
from .model import MCMC
from .tools import get_orientation


if __name__ == "__main__":
    # Example
    path_layout = "/Users/hlibisev/Desktop/hic/chr1.layout.txt"
    path_lens = "/Users/hlibisev/Desktop/hic/comp18_lens.tsv"
    path_pairs = "/Users/hlibisev/Desktop/hic/pairs18.txt"

    # longest_contig
    pairs, contigs, id_contig, longest_contig = get_contigs_and_pairs(path_layout, path_lens, path_pairs, long_contig=True)

    print("Estimation of density...")
    P, f = density(longest_contig)
    print("Estimation of density is done")

    print("MCMC is running...")
    get_orientation([1 for i in range(len(contigs))], pairs, contigs)
    MCMC(pairs, contigs, P, 80, True, n_chains=1)
    print("Have found follow orientation:", [contigs[i].o for i in range(len(contigs))])

    with open("../data/final.layout.txt", "w") as file:
        sign = lambda x: "+" if x == 1 else "-"
        file.write(",".join([contig.name + sign(contig.o) for contig in contigs]))
    print("Result has been saved")