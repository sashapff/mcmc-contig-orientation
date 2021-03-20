from load import get_contigs_and_pairs
from prob import density
from model import MCMC
from tools import get_orientation
import numpy as np


if __name__ == "__main__":
    print("Start!")
    chr_ind = sys.argv[1]

    # Example

    path_layout = "/GWSPH/groups/cbi/Users/pavdeyev/HiCProject/layouts/chr" + chr_ind + "1.layout.txt"
    path_lens = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/contig.length.txt"
    path_pairs = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/pairs/chr_pairs" + chr_ind + ".txt"

    # longest_contig
    pairs, contigs, id_contig, longest_contig = get_contigs_and_pairs(path_layout, path_lens, path_pairs, long_contig=True)

    correct_contigs = [contig.o for contig in contigs]

    print("Estimation of density...")
    P, f = density(longest_contig)
    print("Estimation of density is done")

    print("MCMC is running...")
    get_orientation([1 for i in range(len(contigs))], pairs, contigs)
    MCMC(pairs, contigs, P, 100, n_chains=1)
    print("Have found follow orientation:", [contigs[i].o for i in range(len(contigs))])


    with open("/lustre/groups/cbi/Users/aeliseev/aivanova/data/final" + chr_ind + ".layout.txt", "w") as file:
        sign = lambda x: "+" if x == 1 else "-"
        file.write(",".join([contig.name + sign(contig.o) for contig in contigs]))
    print("Result has been saved")

    print("Count correctness...")
    correct_number = np.array([contig.o == correct_contigs[id_contig[contig.name]] for contig in contigs]).sum()
    print(f"{correct_number}/{len(contigs)} contigs were oriented correctly ({correct_number / len(contigs) * 100}%)")