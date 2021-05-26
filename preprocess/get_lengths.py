import h5py

if __name__ == "__main__":
    chr_indexes = ['MT', 'X']
    for i in range(1, 23):
        chr_indexes.append(str(i))
    contigs2chr = {}
    output_files = {}

    path_layouts = "/GWSPH/groups/cbi/Users/pavdeyev/HiCProject/layouts"
    path_output = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/contig_length"

    matrix_filename = "/GWSPH/groups/cbi/Users/pavdeyev/HiCProject/chm13.draft_v0.9.matrix.cool"

    print('Reading layouts...')
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
        output_files[ind] = open(path_output + "/contig.length." + ind + ".txt", "w")

    print(len(contigs2chr))

    print('Reading matrix...')
    with h5py.File(matrix_filename, "r") as f:
        for i in range(len(f['chroms']['length'][:])):
            contig_name = (f['chroms']['name'][:][i]).decode("utf-8")
            string_to_write = contig_name + '\t' + str(f['chroms']['length'][:][i]) + '\n'
            if contig_name in contigs2chr:
                output_files[contigs2chr[contig_name]].write(string_to_write)
            else:
                print(contig_name)

    for ind in chr_indexes:
        output_files[ind].close()
