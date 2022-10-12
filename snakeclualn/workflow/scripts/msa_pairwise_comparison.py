""" Script pour calculer le pourcentage d'identité et de 
    similiditude (àpartir d'une matrice: BLOSUM 62 par exemple)
    à partir d'un alignement multiple au format fasta
----------------------------------------------------------------
Stagiaire: BARRY Fatoumata Binta M1 BI Université Paris Diderot
Tutrice : DUPRAT Elodie 
Année : 2019-2020 
----------------------------------------------------------------

"""
# Modules utiles

from concurrent.futures import process
import sys
import os
import numpy as np
import pandas as pd
import argparse
import time
import multiprocessing 

from tqdm import tqdm       


    

import log

logger = log.setlogger("pairwise_comparison")

# Fonctions
def parse_dict(old_dict):
    new_dict={}
    for ks,i in old_dict.items():
        k1,k2=ks
        if not k1 in new_dict.keys():
            new_dict[k1]={}
        if not k2 in new_dict.keys():
            new_dict[k2]={}
        new_dict[k1].update({k2:i})
        new_dict[k2].update({k1:i})
        new_dict[k1].update({k1:100})
        new_dict[k2].update({k2:100})
    return new_dict

def parse_dict_distance(old_dict):
    new_dict={}
    for ks,i in old_dict.items():
        k1,k2=ks
        if not k1 in new_dict.keys():
            new_dict[k1]={}
        if not k2 in new_dict.keys():
            new_dict[k2]={}
        new_dict[k1].update({k2:i})
        new_dict[k2].update({k1:i})
        new_dict[k1].update({k1:0})
        new_dict[k2].update({k2:0})
    return new_dict

def parse_scoring_matrix(matrice):
    """ 
        load matrix and parse it into lists

    Parameters
    ----------
    matrix : txt file
        Scoring matrix

    Returns
    -------
    lists
        similitude: list of similar aa (score > 0)
        diff: list of aa with a score egal to 0
    
    """
    # global SIMILAR_AA
    with open(matrice, "r") as mat:
        file = []
        for line in mat:
            if "#" in line:
                continue
            else:
                file.append(line[:-1].split())

    dico_mat = {}
    for i in range(1, len(file)):
        for j in range(len(file[0])):
            aa = (file[0][j], file[i][0])
            dico_mat[aa] = file[i][j+1]

    SIMILAR_AA = []
    diff = []
    for i in dico_mat:
        if int(dico_mat[i]) > 0:
            SIMILAR_AA.append(i)
        else:
            diff.append(i)
    
    return SIMILAR_AA, diff


def parse_fasta(fasta_fichier):
    """ Lit un fichier au format fasta.

    Parametre
    ----------
    fasta_fichier : argument (sys.arg[1])
        file_fasta du fichier au format fasta

    Returns
    -------
    dict
        Dictionnaire avec comme clé l'identifiant de la proteine
        et comme valeur sa sequence proteique sous forme de
        chaine de caractère.

    """
    prot_dict = {}
    with open(fasta_fichier, "r") as fasta:
        prot_id = ""
        for line in fasta.readlines():   
            line=line.strip()     
            if line.startswith(">"):            
                prot_id = line.split()[0].replace(">","")                                
                prot_dict[prot_id] = ""
            else:
                prot_dict[prot_id] += line.strip()
        return(prot_dict)

def match_missmatch_count(dico_fasta):
    """
        return a tuple (id1,id2,count,len1,len2)
    """

def _pairwise_comparison(seq1,seq2,similar_aa):
    """
        match , identity match and gap counts
    """

    match=0
    gap=0
    ident=0
    similitude=0

    # retrieve sequence length without gap
    len_seq1 = len(seq1.replace("-",""))
    len_seq2 = len(seq2.replace("-",""))

    for i,j in zip(seq1,seq2):
        if i!="-" and j!="-":
            match+=1 #match    
            if i==j:
                ident+=1                
        elif i == "-" or j == "-":
            gap += 1 #gap             
        else:
            sys.exit()
        if (i,j) in similar_aa or (j,i) in similar_aa:
            similitude += 1

    shortest_seq_length = min([len_seq1,len_seq2])
    
    try:
        perc_of_id  = round(100*ident / shortest_seq_length,1)
        identity_distance = round(100-100*ident / shortest_seq_length, 1)/100
        perc_of_sim = round(100*similitude / shortest_seq_length, 1)
        similarity_distance = round(100-100*similitude / shortest_seq_length, 1)/100
        # compute coverage (number of aa aligned between seq1 and seq2 even if it's a missmatch)                                                
        coverage_1_2 = match/len_seq1
        coverage_2_1 = match/len_seq2           
    except ZeroDivisionError:        
        return None, None , None ,  None , None , None , None , None , None , None , None        
    return len_seq1, len_seq2 , len(seq2),  match , perc_of_id , identity_distance , perc_of_sim , similarity_distance , coverage_1_2 , coverage_2_1 , gap


def pairwise_comparison(args):
    
    set1,set2,similar_aa,combi,pos=args

    data=[]
    
    logger.info("Running batch {} vs batch {}. ".format(combi[0],combi[1]))
    
    for id1,seq1 in tqdm( set1.items() , position=pos, leave=True , desc = "Batches : {}".format(str(combi))) :        
    
        for id2,seq2 in set2.items():
            # _pairwise_comparison (compare each position between seq1 and seq2) - return match count, identity count and gap count as tuple                      
            d = (id1,id2)   
               
            metrics = _pairwise_comparison(seq1,seq2,similar_aa)                
            d += metrics      
            data.append(d)
                        
        set2.pop(id1, None)

    return data


def split_dict(d, n):
    keys = list(d.keys())
    batches = []
    for i in range(0, len(keys), n):
        batches.append( {k: d[k] for k in keys[i: i + n]} )
    return batches


def run_comparison( dico_fasta , similar_aa):
    batch_size = 100

    dico_fasta = dict(sorted(dico_fasta.items()))

    batches = split_dict(dico_fasta,batch_size)

    logger.info("{} comparisons ( N * (N-1)  where N={}).".format(
            len(dico_fasta) * (len(dico_fasta) - 1 ),
            len(dico_fasta)
        )
    )

    st = time.time()
    data=[]
    processes = []
    combi = []
    pos = 0
    for b1 in range(0,len(batches),1):
        for b2 in range(0,len(batches),1):            
            if (b1,b2) not in combi or (b2,b1) not in combi:
                combi+= [(b1,b2),(b2,b1)]
                logger.info("batches {}-{}".format(b1,b2))
                processes.append(
                    ( batches[b1].copy() , batches[b2].copy() , similar_aa , (b1,b2) , pos)
                )
            pos += 1
    
    threads = multiprocessing.cpu_count() - 1               
    logger.info("Compare sequences by batch [n={} seq] with {} threads".format(batch_size,threads))
    with multiprocessing.Pool(threads) as p:  
              
        data += p.imap(pairwise_comparison , processes)

    et = time.time()
    elapsed_time = et - st
    print('Execution time in multiprocessing mode :', elapsed_time, 'seconds')
    # st = time.time()
    # seq1_d = dico_fasta.copy()
    # seq2_d = dico_fasta.copy()
    # data=[]
    # for id1,seq1 in seq1_d.items():        
    #     for id2,seq2 in seq2_d.items():
    #         # _pairwise_comparison (compare each position between seq1 and seq2) - return match count, identity count and gap count as tuple                      
    #         d = (id1,id2)
    #         try:
    #             metrics = _pairwise_comparison(seq1,seq2,similar_aa)
    #             d += metrics                                            
    #             data.append(d)
    #         except ZeroDivisionError:
    #             pass
                
    #     seq2_d.pop(id1, None)
    # et = time.time()
    # elapsed_time = et - st
    # print('Execution time in sequential mode :', elapsed_time, 'seconds')
    # data is a list of tuple with the following value for each comparison : 
    #   "seq1 identifier"
    #   "seq2 identifier"
    #   "seq1 length"
    #   "seq2 length"
    #   "aln length"
    #   "seq1 coverage"
    #   "seq2 coverage"
    #   "seq1 identifier"
    #   "percentage of identity" where pid here is the number of correct match divide by the length of the shortest sequence.
    #   "number of gap"    
    return data[0]



def make_matrix(df,metric):
    pairwise_dict = {}
    for i,j in df.iterrows():            
        if j["id1"] not in pairwise_dict:
            pairwise_dict[j["id1"]] = {}
        if j["id2"] not in pairwise_dict:
            pairwise_dict[j["id2"]] = {}
        pairwise_dict[j["id1"]].update({j["id2"] : j[metric]})        
        pairwise_dict[j["id2"]].update({j["id1"] : j[metric]}) 
    return pd.DataFrame.from_dict(pairwise_dict)


def parseargs():
    parser = argparse.ArgumentParser(
        prog=''
    )
    #parser.add_argument('-',nargs='+',default=sys.stdin)
    parser.add_argument("-i",'--infile',help="fasta file")
    parser.add_argument('-m','--matrice',required=True,help="scoring matrice e.g BLOSUM62")
    parser.add_argument('--distance',action="store_true",help="save also distance matrices from identity and similarity scores")
    parser.add_argument('-o',default=sys.stdout,help="output directory")
    parser.add_argument('--log',default=sys.stdout,help="output directory")
    args = parser.parse_args()
    return args

def parsesnake():
    args = argparse.Namespace(infile=str(snakemake.input.msa),
            matrice = str(snakemake.params.matrix),
            distance = str(snakemake.params.distance),
            o = str(snakemake.params.outdir),
            log = str(snakemake.log)
        )
    return args

def generate_outfilename(outdir,filename):
    if outdir is sys.stdout:
        return sys.stdout
    else:
        return os.path.join(outdir,filename + ".tsv")

if __name__ == '__main__':   
    if 'snakemake' in globals():
        args = parsesnake()
    else:
        args = parseargs()

    if args.log != sys.stdout:
        logdir = os.path.dirname(os.path.abspath(args.log))            
        
        os.makedirs(logdir,exist_ok=True)
        
        logger.addHandler(
            log.file_handler(args.log,"INFO")
        )
    
    else:
        logger.addHandler(
            log.stream_handler("INFO")
        )

    file_fasta = args.infile
    matrix = args.matrice
    outdir = args.o

    if isinstance(outdir,str):
        os.makedirs(outdir,exist_ok=True)


    # load matix
    similar_aa , _ = parse_scoring_matrix(matrix)

    # fasta into in dictionnary
    dico_fasta = parse_fasta(file_fasta)

    # si fasta > 1 do something
    if len(dico_fasta) > 1:
        logger.info("Fasta file is suitable for pairwise comparison (n>1) ")
        
        # 1 : pairwise_iteration
        pairwise_df_flat = pd.DataFrame(
            run_comparison( dico_fasta , similar_aa )
            )
        
        #len_seq1, len_seq2 , len(seq2),  match , perc_of_id , identity_distance , perc_of_sim , similarity_distance , coverage_1_2 , coverage_2_1 , gap
        pairwise_df_flat.columns="id1,id2,len_seq1,len_seq2,len_aln,match,identity_percent,identity_distance,similarity_percent,similarity_distance,coverage_1_2,coverage_2_1,gap_count".split(",")
        
        logger.info("Building matrices.")
        ipdf = make_matrix(pairwise_df_flat,"identity_percent")
        iddf = make_matrix(pairwise_df_flat,"identity_distance")
        spdf = make_matrix(pairwise_df_flat,"similarity_percent")
        sddf = make_matrix(pairwise_df_flat,"similarity_distance")
        logger.info("Done.")

        ipdf.to_csv( generate_outfilename(outdir, "identity_percent")  ,sep="\t",header=True,index=True)        
        iddf.to_csv( generate_outfilename(outdir, "identity_distance"), sep="\t",header=True,index=True)        
        spdf.to_csv(generate_outfilename(outdir, "similarity_percent"), sep="\t",header=True,index=True)        
        sddf.to_csv( generate_outfilename(outdir, "similarity_distance"), sep="\t",header=True,index=True)                    
        pairwise_df_flat.to_csv( generate_outfilename(outdir,"summary") ,sep="\t",header=True,index=False)
        
    else:
        logger.info("Fasta file is not suitable for pairwise comparison (n==1)")        
        if 'snakemake' in globals():
            logger.info("Creating dummy output file for snakemake")
            open(os.path.join(outdir,"identity_percent.tsv"),'w').close()
            open(os.path.join(outdir,"identity_distance.tsv"),'w').close()
            open(os.path.join(outdir,"similarity_percent.tsv"),'w').close()
            open(os.path.join(outdir,"similarity_distance.tsv"),'w').close()
            open(os.path.join(outdir,"summary.tsv"),'w').close()


    logger.info("End.")

    #     with open(os.path.join())
    #     for i in UPGMA(li_sim, ide):
    # #     #dico_dis_simi = parse_dict(distance)
    # # pd.DataFrame.from_dict(,orient="index").to_csv(os.path.join(outdir,"simi_distance.tsv"),sep="\t",header=True,index=True)
    # # #pd.DataFrame.from_dict(UPGMA(li_iden, ide),orient="index").to_csv(os.path.join(outdir,"ident_distance.tsv"),sep="\t",header=True,index=True)
    