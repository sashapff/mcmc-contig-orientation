import h5py
import numpy as np

if __name__ == "__main__":
    matrix_filename = "/GWSPH/groups/cbi/Users/pavdeyev/HiCProject/chm13.draft_v0.9.matrix.cool"
    output_file = open("/lustre/groups/cbi/Users/aeliseev/aivanova/data/contig.length.txt", "w")

    with h5py.File(matrix_filename, "r") as f:
        print(f['resolutions'].keys())
        # for i in range(len(f['resolutions']['50000']['chroms']['length'][:])):
        #     string_to_write = (f['resolutions']['50000']['chroms']['name'][:][i]).decode("utf-8")  + '\t' \
        #                       + str(f['resolutions']['50000']['chroms']['length'][:][i]) + '\n'
        #     output_file.write(string_to_write)

    output_file.close()
