import sys
import argparse

import pandas as pd
from Bio import SeqIO

# try:
#     import snakemake
#     print("get_representatives.py in snakemake mode")
# except ModuleNotFoundError:
#     pass

if __name__ == '__main__':
    try:
        fasta = str(snakemake.input.fasta)
        clufile = str(snakemake.input.clufile)
        outfile = str(snakemake.output)
    except:
        parser = argparse.ArgumentParser(
            prog='get_representatives.py',
            description='extract representative sequences from fasta and a cluster file'
        )
        parser.add_argument("-f",'--fasta'  , type=str,required=True)
        parser.add_argument("-c",'--clufile', type=str,required=True)
        parser.add_argument('-o','--outfile', default=sys.stdout,help="outfile")
        args = parser.parse_args()
        fasta = args.fasta
        clufile = args.clufile
        outfile = args.outfile

    representatives_list = list(set(pd.read_csv(clufile,sep="\t",header=None,index_col=0).index))

    record_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    rep_record_list = []
    for rid,record in record_dict.items():
        if rid in representatives_list:
            rep_record_list.append(record)

    if outfile is not sys.stdout:
        SeqIO.write(rep_record_list , outfile , "fasta")
    else:
        for record in rep_record_list:
            r = ">{} {}\n{}\n".format(record.id , record.description , record.seq )
            print(r)

