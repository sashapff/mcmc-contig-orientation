from tools.load import check_reads, get_contigs_and_pairs

if __name__ == "__main__":
    chr_ind = 3
    path_pairs = f'/lustre/groups/cbi/Users/aeliseev/aivanova/data/pairs/chr_pairs{3}.txt'
    output_path = f'/lustre/groups/cbi/Users/aeliseev/aivanova/data/check_chr{3}.txt'
    # path_pairs = '/Users/alexandra/bioinf/mcmc/chr3/chr_pairs3.txt'
    check_reads(path_pairs, output_path, chr_ind)

    # pairs, contigs, id_contig, longest_contig_pairs, longest_contig = get_contigs_and_pairs(path_layout, path_lens,
    #                                                                                         path_pairs,
    #                                                                                         long_contig=True,
    #                                                                                         min_len=min_contig_length)