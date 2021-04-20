from tools.load import check_reads

if __name__ == "__main__":
    path_pairs = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/pairs/chr_pairs3.txt"
    print(check_reads(path_pairs))