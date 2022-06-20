import os
import sys
import re

files = snakemake.params.files
msas = snakemake.output.msas
tbl = snakemake.output.tbl


# make input_msa
with open(str(msas),'w') as output_handle:
    tbl_dict = {}
    nseq=0
    for i in files:
        filename = os.path.basename(i)
        lseq = " "
        with open(i,'r') as file_handle:
            for l in file_handle.readlines():
                # add sequence to msas output (concat all msa)
                output_handle.write(l.strip()+"\n")
                # count the number of sequence
                if l.startswith('>'):
                    nseq += 1
                    lseq += " {}".format(str(nseq))
        tbl_dict[filename] = lseq

with open(str(tbl),'w') as output_handle:
    for i,j in tbl_dict.items():
        # write clu comp : 1 2 # cluster_1.fasta   
        l = "{}  # {}\n".format(j,i)
        output_handle.write(l)
