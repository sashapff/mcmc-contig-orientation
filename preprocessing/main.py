if __name__ == "__main__":
    chr_indexes = ['MT', 'X']
    for i in range(1, 23):
        chr_indexes.append(str(i))
    contigs2chr = {}
    chr_pairs = {}

    path_layouts = "/GWSPH/groups/cbi/Users/pavdeyev/HiCProject/layouts"
    path_output = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/pairs"

    for ind in chr_indexes:
        if ind != "9":
            filename = path_layouts + "/chr" + ind + ".layout.txt"
        else:
            filename = path_layouts + "/chr9.partial_layout.txt"
        with open(filename, "r") as f:
            lines = f.read().splitlines()
        for line in lines:
            chr_name, contigs = line.split(" ")
            contigs = contigs.split(",")
            for contig in contigs:
                contigs2chr[contig[:-1]] = ind
        chr_pairs[ind] = open(path_output + "/chr_pairs" + ind + ".txt", "w")

    print("Reading pairs...")

    path_pairs = "/GWSPH/groups/cbi/Users/pavdeyev/HiCProject/chm13.draft_v0.9.simplified.nodes.pairs"
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
