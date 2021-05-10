from tools.load import get_contigs_and_pairs
from tools.prob import density
import matplotlib.pyplot as plt

if __name__ == "__main__":
    path_to_output = f'/lustre/groups/cbi/Users/aeliseev/aivanova/data/density'
    # path_to_output = f'/Users/alexandra/bioinf/mcmc/output/{min_contig_length_name}'

    chr_ind = '4'
    path_layout = "/GWSPH/groups/cbi/Users/pavdeyev/HiCProject/layouts/chr" + chr_ind + ".layout.txt"
    path_lens = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/contig_length/contig.length." + chr_ind + ".txt"
    path_pairs = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/pairs/chr_pairs" + chr_ind + ".txt"
    # longest_contig
    pairs, contigs, id_contig, longest_contig_pairs, longest_contig = get_contigs_and_pairs(path_layout, path_lens,
                                                                                            path_pairs,
                                                                                            long_contig=True,
                                                                                            min_len=0,
                                                                                            from_one_contig=True)

    P, f = density(longest_contig_pairs)

    p_range = range(110000)
    plt.title('The density of distance between pair of reads')
    plt.plot(p_range, [P(x) for x in p_range])
    plt.xlabel('x, distance')
    plt.ylabel('p(x), density in distance x')
    plt.legend()
    plt.savefig(f'{path_to_output}/density_100k.png')

    plt.clf()

    p_range = range(5000)
    plt.title('The density of distance between pair of reads')
    plt.plot(p_range, [P(x) for x in p_range])
    plt.xlabel('x, distance')
    plt.ylabel('p(x), density in distance x')
    plt.legend()
    plt.savefig(f'{path_to_output}/density_5k.png')

    plt.clf()

    p_range = range(5000)
    plt.title('The density of distance between pair of reads')
    plt.plot(p_range, [f(x) for x in p_range])
    plt.xlabel('x, distance')
    plt.ylabel('p(x), density in distance x')
    plt.legend()
    plt.savefig(f'{path_to_output}/density_f.png')
