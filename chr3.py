from tools.load import check_reads

if __name__ == "__main__":
    chr_ind = 3
    path_pairs = f'/lustre/groups/cbi/Users/aeliseev/aivanova/data/pairs/chr_pairs{3}.txt'
    output_path = f'/lustre/groups/cbi/Users/aeliseev/aivanova/data/check_chr{3}.txt'
    # path_pairs = '/Users/alexandra/bioinf/mcmc/chr3/chr_pairs3.txt'
    check_reads(path_pairs, output_path, chr_ind)