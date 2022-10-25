#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import os
import sys
import numpy as np
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse
import re

if __name__ == '__main__':
    try:
        fasta = str(snakemake.input)
        outfile = str(snakemake.output)
    except:
        parser = argparse.ArgumentParser(
            prog='get_seq.py',
            description='get sequence length from a fasta file '
        )
        #parser.add_argument('-',nargs='+',default=sys.stdin)
        parser.add_argument("-i",'--infile', type=str,required=True)
        parser.add_argument('-o',default=sys.stdout,help="outfile")
        args = parser.parse_args()
        fasta = args.infile
        outfile = args.o
    

    data=[]
    with open(fasta,'r') as stream:
        for r in SeqIO.parse(stream, 'fasta'):    
            data.append((r.id, len(str(r.seq)) , str(r.seq) ))


    if outfile is sys.stdout:
        for i in data:
            print(i[0].strip(),"\t",i[1] , i[2])
    else:
        with open(str(outfile),"w") as fout:
            for i in data:
                
                fout.write("{}\t{}\t{}".format(i[0].strip(),i[1] , i[2] ))
                # else:
                    # fout.write("{}\t{}".format(i[0].strip(),i[1]))

            
            
            
            
            
            
            
            
            
            

