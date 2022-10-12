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

import sys
import os
import numpy as np
import pandas as pd
import argparse
try:
    from tqdm import tqdm
    print("tqdm loaded")
    tqdm([i for i in range(0,100000)])
except ModuleNotFoundError:
    print("error while loading tqdm")
    
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

    similitude = []
    diff = []
    for i in dico_mat:
        if int(dico_mat[i]) > 0:
            similitude.append(i)
        else:
            diff.append(i)
    return similitude, diff


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

    perc_of_id  = round(100*ident/(len(seq2)-gap),1)
    identity_distance = round(100-100*ident/(len(seq2)-gap), 1)/100
    perc_of_sim = round(100*similitude/(len(seq2)-gap), 1)
    similarity_distance = round(100-100*similitude/(len(seq2)-gap), 1)/100
    # compute coverage (number of aa aligned between seq1 and seq2 even if it's a missmatch)                                            
    try:
        coverage_1_2 = match/len_seq1
        coverage_2_1 = match/len_seq2            
    except ZeroDivisionError:
        coverage_1_2 = None
        coverage_2_1 = None

    return len_seq1, len_seq2 , len(seq2),  match , perc_of_id , identity_distance , perc_of_sim , similarity_distance , coverage_1_2 , coverage_2_1 , gap

def pairwise_comparison(dico_fasta,matrix):
    dico_fasta = dict(sorted(dico_fasta.items()))
    similar_aa, _  = parse_scoring_matrix(matrix)
    seq1_d = dico_fasta.copy()
    seq2_d = dico_fasta.copy()
    data=[]
    for id1,seq1 in tqdm(seq1_d.items()):        
        for id2,seq2 in seq2_d.items():
            # _pairwise_comparison (compare each position between seq1 and seq2) - return match count, identity count and gap count as tuple            
            
            d = (id1,id2)
            metrics = _pairwise_comparison(seq1,seq2,similar_aa)
            d += metrics
            
            try:                
                data.append(d)                
            except ZeroDivisionError:                  
                data.append(d)
                
        seq2_d.pop(id1, None)
    # data is a list of tuple with the following value for each comparison : 
    #   "seq1 identifier"
    #   "seq2 identifier"
    #   "seq1 length"
    #   "seq2 length"
    #   "aln length"
    #   "seq1 coverage"
    #   "seq2 coverage"
    #   "seq1 identifier"
    #   "percentage of identity" where pid here is the number of correct match divide by the length of the alignment minus gap.
    #   "number of gap"
    return data



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
    args = parser.parse_args()
    return args

def parsesnake():
    args = argparse.Namespace(infile=str(snakemake.input.msa),
            matrice = str(snakemake.params.matrix),
            distance = str(snakemake.params.distance),
            o = str(snakemake.params.outdir)
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

    file_fasta = args.infile
    matrice = args.matrice
    outdir = args.o

    if isinstance(outdir,str):
        os.makedirs(outdir,exist_ok=True)


    # fasta into in dictionnary
    dico_fasta = parse_fasta(file_fasta)
    
    # si fasta > 1 do something
    if len(dico_fasta) > 1:
        # 1 : pairwise_iteration
        pairwise_df_flat = pd.DataFrame(pairwise_comparison(dico_fasta , matrice ))
        
        #len_seq1, len_seq2 , len(seq2),  match , perc_of_id , identity_distance , perc_of_sim , similarity_distance , coverage_1_2 , coverage_2_1 , gap
        pairwise_df_flat.columns="id1,id2,len_seq1,len_seq2,len_aln,match,identity_percent,identity_distance,similarity_percent,similarity_distance,coverage_1_2,coverage_2_1,gap_count".split(",")
        

        ipdf = make_matrix(pairwise_df_flat,"identity_percent")
        iddf = make_matrix(pairwise_df_flat,"identity_distance")
        spdf = make_matrix(pairwise_df_flat,"similarity_percent")
        sddf = make_matrix(pairwise_df_flat,"similarity_distance")

    
        print("Identity (%): \n")
        ipdf.to_csv( generate_outfilename(outdir, "identity_percent")  ,sep="\t",header=True,index=True)
        print("Identity (distance): \n")
        iddf.to_csv( generate_outfilename(outdir, "identity_distance"), sep="\t",header=True,index=True)
        print("\nSimilarity (%) : \n")
        spdf.to_csv(generate_outfilename(outdir, "similarity_percent"), sep="\t",header=True,index=True)
        print("\nSimilarity (distance) : \n")
        sddf.to_csv( generate_outfilename(outdir, "similarity_distance"), sep="\t",header=True,index=True)
                    
        pairwise_df_flat.to_csv( generate_outfilename(outdir,"summary") ,sep="\t",header=True,index=False)

    else:
        open(os.path.join(outdir,"identity_percent.tsv"),'w').close()
        open(os.path.join(outdir,"identity_distance.tsv"),'w').close()
        open(os.path.join(outdir,"similarity_percent.tsv"),'w').close()
        open(os.path.join(outdir,"similarity_distance.tsv"),'w').close()
        open(os.path.join(outdir,"summary.tsv"),'w').close()




    #     with open(os.path.join())
    #     for i in UPGMA(li_sim, ide):
    # #     #dico_dis_simi = parse_dict(distance)
    # # pd.DataFrame.from_dict(,orient="index").to_csv(os.path.join(outdir,"simi_distance.tsv"),sep="\t",header=True,index=True)
    # # #pd.DataFrame.from_dict(UPGMA(li_iden, ide),orient="index").to_csv(os.path.join(outdir,"ident_distance.tsv"),sep="\t",header=True,index=True)
    