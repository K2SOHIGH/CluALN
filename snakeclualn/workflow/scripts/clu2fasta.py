#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import os
import sys
import numpy as np
import re
import gzip 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse

def fasta2records(fasta):
    try:
        records = SeqIO.to_dict(SeqIO.parse(open(fasta,"rt"), 'fasta'))    
        return records # { k.split("|")[1] : v if re.match('lcl|ref',k,re.IGNORECASE) else k for k,v in records.items()}
    except ValueError: #duplicated keys
        records_dict = dict()
        for record in SeqIO.parse(fasta, 'fasta'):
            if record.id not in records_dict.keys():
                #new record
                records_dict[record.id]=record
            else:
                #duplicate entry
                #check sequence
                if record.seq == records_dict[record.id].seq:
                    #same sequence
                    continue
                else:
                    raise ValueError("{} have mutliple entry but sequences are different ...".format(record.id))
        return records_dict

        
    #return { re.sub("\|.*", "" ,k.replace("lcl|","")) : v  for k,v in records.items()}

def write_fasta(rep,l_records,outfile):
    if outfile is not sys.stdout:
        os.makedirs(outfile,exist_ok=True)
        outfile = os.path.join(outfile,rep+".fasta")
        with open(str(outfile), "w") as handle:
            SeqIO.write(l_records, handle, "fasta")
    else:
        sys.stdout.write("Cluster rep : {}\n".format(rep))          
        for r in l_records:  
            sys.stdout.write(">{}\t{}\n".format(r.id,r.description))
            sys.stdout.write("{}\n".format(str(r.seq)))

def parseargs():
    parser = argparse.ArgumentParser(
        prog='filterdomtbl.py ',
        description='extract hits bound from hmmsearch output (no filtering)'
    )

    #"python3 clu2fasta.py -i {input.clu} -f {input.fa} -o {params.outdir} && touch {output}"
    parser.add_argument('-i','--input',dest="clu",required=True,help="cluster tabular file from mmseqs clustering")
    parser.add_argument('-f','--fasta',dest="fasta",required=True,help="fasta file used for clustering")
    parser.add_argument('-o','--out',default=sys.stdout,help="outdir - where clu_<identifier>.fasta will be stored   [sys.stdout]")
    args = parser.parse_args()
    return args

def parsesnake():
    args = argparse.Namespace(clu=str(snakemake.input.clu),
        fasta = str(snakemake.input.fa),
        out = str(snakemake.params.outdir)
        )
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        args = parsesnake()
    else:
        args = parseargs()


    df = pd.read_csv(args.clu,sep="\t",header=None)
    df.columns=["rep","seq"]

    # get all sequence into a records dict
    records=fasta2records(args.fasta)
    # sequence identifier
    records_id_l = records.keys()
    cluster_id = 1
    #groupby cluster and subset fasta to cluster_fasta
    for rep, sdf in df.groupby("rep"):
        members = list(sdf.seq)
        clu_rec=[]
        for i in members:
            if i in records:
                clu_rec.append(records[i])
        write_fasta("cluster_"+str(cluster_id),clu_rec,args.out)
        cluster_id+=1
    if 'snakemake' in globals():
        open(os.path.join(args.out,"tofastaa.done"),'w').close()