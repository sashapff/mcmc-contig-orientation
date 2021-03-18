if __name__ == "__main__":
    chr_indexes = [1]
    contigs2chr = {}
    chr_pairs = {}

    for ind in chr_indexes:
        with open("../data/chr" + str(ind) + ".layout.txt", "r") as f:
            lines = f.read().splitlines()
        for line in lines:
            chr_name, contigs = line.split(" ")
            contigs = contigs.split(",")
            for contig in contigs:
                contigs2chr[contig[:-1]] = ind
        chr_pairs[ind] = open("../data/chr_pairs" + str(ind) + ".txt", "w")

    path_pairs = "../data/pairs18.txt"
    with open(path_pairs, "r") as f:
        line = f.readline()
        for line in f:
            s = line.strip().split("\t")
            contig_1 = s[1]
            contig_2 = s[3]
            if contig_1 in contigs2chr and contig_2 in contigs2chr and contigs2chr[contig_1] == contigs2chr[contig_2]:
                chr_pairs[contigs2chr[contig_1]].write(line)

    for ind in chr_indexes:
        chr_pairs[ind].close()