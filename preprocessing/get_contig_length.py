import h5py
import numpy as np

if __name__ == "__main__":
    matrix_filename = "/GWSPH/groups/cbi/Users/pavdeyev/HiCProject/chm13.draft_v0.9.matrix.cool"
    output_file = open("/lustre/groups/cbi/Users/aeliseev/aivanova/data/contig.length.txt", "w")

    # path = "/Users/alexandra/bioinf/mcmc/data/"
    #
    # matrix_filename = path + "mat18_50k.cool"
    # output_file = open(path + "contig_length.txt", "w")

    with h5py.File(matrix_filename, "r") as f:
        for i in range(len(f['chroms']['length'][:])):
            string_to_write = (f['chroms']['name'][:][i]).decode("utf-8")  + '\t' \
                              + str(f['resolutions']['50000']['bins']['length'][:][i]) + '\n'
            output_file.write(string_to_write)

    output_file.close()
