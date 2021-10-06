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

import sys, os
import numpy as np
import pandas as pd
import argparse


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

def mat_simi(matrice):
    """ Lit une matrice de score.

    Parametre
    ----------
    matrice : fichier txt
        Matrice de score

    Returns
    -------
    listes
        similitude: acides aminés similaires dont le score > 0
        diff: acides aminés non similaires
    
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


def lit_fasta(fasta_fichier):
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

def iter_seq(seq1,seq2):
    cpt=0
    gap=0
    ident=0
    for i,j in zip(seq1,seq2):
        if i!="-" and j!="-":
            cpt+=1 #match    
            if i==j:
                ident+=1                
        elif i == "-" or j == "-":
            gap += 1 #gap             
        else:
            sys.exit()
    return cpt,ident,gap

#dico_ident[ide[i],ide[j]] = round(100*ident/(l-gap),1)

def pairwise_iteration(dico_fasta):
    seq1_d = dico_fasta.copy()
    seq2_d = dico_fasta.copy()
    data=[]
    for id1,seq1 in seq1_d.items():
        seq2_d.pop(id1, None)
        for id2,seq2 in seq2_d.items():
            mm_count,ident,gap = iter_seq(seq1,seq2)
            len_seq1 = len(seq1.replace("-",""))
            len_seq2 = len(seq2.replace("-",""))
            try:
                coverage_1_2 = mm_count/len_seq1
                coverage_2_1 = mm_count/len_seq2
                data.append((id1,id2,len_seq1,len_seq2,len(seq2),coverage_1_2,coverage_2_1,round(100*ident/(len(seq2)-gap),1),gap))
            except ZeroDivisionError:
                coverage_1_2 = None
                coverage_2_1 = None
                data.append((id1,id2,len_seq1,len_seq2,len(seq2),coverage_1_2,coverage_2_1,None,gap))
    return data


def pourcentage_identite(n,l,sequences,ide):
    """ Calcul le pourcentage d'identité par paire de
        séquences alignées

    Parametres
    ----------
    n : integer
        file_fastabres de sequences alignées
    l : integer
        longueur de l'alignement
    sequences: list
        liste des sequences alignées
    ide: list
        liste des identifiants des protéines
    Returns
    -------
    dict
        Dictionnaire avec comme clé les identifiants par paire
        et comme valeur son pourcentage d'identité
    
    """
    dico_ident = {}
    distan_iden={}
    for i in range(n-1):
        for j in range(i+1,n):
            #paire séquences i,j
            ident = gap = 0
            
            for c in range(l):
                #gap dans l'alignement
                if sequences[i][c] == "-" or sequences[j][c] == "-":
                    gap += 1
                # identite dans l'alignement
                elif sequences[i][c] == sequences[j][c]:
                    ident += 1
            try:
                distan_iden[ide[i],ide[j]] = round(100-100*ident/(l-gap), 1)/100
                dico_ident[ide[i],ide[j]] = round(100*ident/(l-gap),1)
            except ZeroDivisionError:
                distan_iden[ide[i],ide[j]] = 0
                dico_ident[ide[i],ide[j]] = 0
            #if dico_ident[ide[i],ide[j]] is None:
            #if "AFY49796.1" in [ide[i],ide[j]] and "PSO93175.1" in [ide[i],ide[j]]:
                #print(dico_ident[ide[i],ide[j]])
    list_dist=[]
    for i in range(n):
        list_dist.append([])
        for j in range(i):
            if (ide[i],ide[j]) in distan_iden.keys():
                list_dist[i].append(distan_iden[ide[i], ide[j]])
            elif (ide[j],ide[i]) in distan_iden.keys():
                list_dist[i].append(distan_iden[ide[j], ide[i]])

    return dico_ident, distan_iden, list_dist


def pourcentage_similitude(n,l,sequences,ide,similitude):
    """ Calcul le pourcentage de similitude par paire de
        séquences alignées

    Parametres
    ----------
    n : integer
        file_fastabres de sequences alignées
    l : integer
        longueur de l'alignement
    sequences: list
        liste des sequences alignées
    ide: list
        liste des identifiants des protéines
    similitude: list
        liste des acides aminés similaires sous forme de tuple ('A', 'A')
    Returns
    -------
    dict
        Dictionnaire avec comme clé les identifiants par paire
        et comme valeur son pourcentage d'identité
    
    """    
    
    dico_simi = {}
    distan_simi={}
    for i in range(n-1):
        for j in range(i+1,n):
            #paire séquences i,j
            simi=dif=gap=0
            for c in range(l):
                if sequences[i][c] == "-" or sequences[j][c] == "-":
                    gap+=1
                # similitude dans l'alignement
                elif (sequences[i][c],sequences[j][c]) in similitude or ((sequences[j][c],sequences[i][c]) in similitude):
                    simi+=1
                else:
                    dif+=1
            #print(l-gap,l,gap)
            try:
                distan_simi[ide[i],ide[j]] = round(100-100*simi/(l-gap), 1)/100
                dico_simi[ide[i],ide[j]] = round(100*simi/(l-gap), 1)
            except ZeroDivisionError:
                distan_simi[ide[i],ide[j]] = 1
                dico_simi[ide[i],ide[j]] = 0
            
    li=[]
    for i in range(n):
        li.append([])
        for j in range(i):
            if (ide[i],ide[j]) in distan_simi.keys():
                li[i].append(distan_simi[ide[i],ide[j]])
            elif (ide[j],ide[i]) in distan_simi.keys():
                li[i].append(distan_simi[ide[j],ide[i]])

    return dico_simi, distan_simi, li



# copyright https://github.com/lex8erna/UPGMApy

""" 			DEBUT			"""

# A Quick Implementation of UPGMA (Unweighted Pair Group Method with Arithmetic Mean)

# lowest_cell:
#   Locates the smallest cell in the table
def lowest_cell(table):
    # Set default to infinity
    min_cell = float("inf")
    x, y = -1, -1

    # Go through every cell, looking for the lowest
    for i in range(len(table)):
        for j in range(len(table[i])):
            if table[i][j] < min_cell:
                min_cell = table[i][j]
                x, y = i, j

    # Return the x, y co-ordinate of cell
    return x, y


# join_labels:
#   Combines two labels in a list of labels
def join_labels(labels, a, b):
    # Swap if the indices are not ordered
    if b < a:
        a, b = b, a

    # Join the labels in the first index
    labels[a] = "(" + labels[a] + "," + labels[b] + ")"

    # Remove the (now redundant) label in the second index
    del labels[b]


# join_table:
#   Joins the entries of a table on the cell (a, b) by averaging their data entries
def join_table(table, a, b):
    # Swap if the indices are not ordered
    if b < a:
        a, b = b, a

    # For the lower index, reconstruct the entire row (A, i), where i < A
    row = []
    for i in range(0, a):
        row.append((table[a][i] + table[b][i])/2)
    table[a] = row
    
    # Then, reconstruct the entire column (i, A), where i > A
    #   Note: Since the matrix is lower triangular, row b only contains values for indices < b
    for i in range(a+1, b):        
        table[i][a] = (table[i][a]+table[b][i])/2
        
    #   We get the rest of the values from row i

    for i in range(b+1, len(table)):
        table[i][a] = (table[i][a]+table[i][b])/2
        # Remove the (now redundant) second index column entry
        del table[i][b]

    # Remove the (now redundant) second index row
    del table[b]


# UPGMA:
#   Runs the UPGMA algorithm on a labelled table
def UPGMA(table, labels):
    # Until all labels have been joined...
    while len(labels) > 1:
        # Locate lowest cell in the table
        x, y = lowest_cell(table)

        # Join the table on the cell co-ordinates
        join_table(table, x, y)

        # Update the labels accordingly
        join_labels(labels, x, y)

    # Return the final label
    return labels[0]

""" 			FIN			"""

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


if __name__ == '__main__':   
    if 'snakemake' in globals():
        args = parsesnake()
    else:
        args = parseargs()

    file_fasta = args.infile
    matrice = args.matrice
    outdir = args.o
    dico_fasta = lit_fasta(file_fasta)
    pairwise_df_flat = pd.DataFrame(pairwise_iteration(dico_fasta))
    
    pairwise_df_flat.columns="id1,id2,len_seq1,len_seq2,len_aln,coverage_1_2,coverage_2_1,identity_percent,gap_count".split(",")

    print(pairwise_df_flat)
    sequences = list(dico_fasta.values())
    ide = list(dico_fasta.keys())
    n = len(sequences)
    l = len(sequences[0])

 
    similitude = mat_simi(matrice)[0]
    (dico_simi, distance, li_sim) = pourcentage_similitude(n,l,sequences,ide,similitude)
    (dico_ident, dist, li_iden) = pourcentage_identite(n,l,sequences,ide)


    if isinstance(outdir,str):
        os.makedirs(outdir,exist_ok=True)


    dico_ident = parse_dict(dico_ident)
    dico_simi = parse_dict(dico_simi)

    
    df = pd.DataFrame.from_dict(dico_ident,orient="index")
    
    if outdir is sys.stdout:
        print("Identity : \n")
        pd.DataFrame.from_dict(dico_ident,orient="index").to_csv(outdir,sep="\t",header=True,index=True)
        print("\nSimilarity : \n")
        pd.DataFrame.from_dict(dico_simi,orient="index").to_csv(outdir,sep="\t",header=True,index=True)
        print("\nSummary : \n")
        pairwise_df_flat.to_csv(outdir,sep="\t",header=True,index=False)
    else:        
        pd.DataFrame.from_dict(dico_ident,orient="index").to_csv(os.path.join(outdir,"identity.tsv"),sep="\t",header=True,index=True)
        pd.DataFrame.from_dict(dico_simi,orient="index").to_csv(os.path.join(outdir,"similarity.tsv"),sep="\t",header=True,index=True)
        pairwise_df_flat.to_csv(os.path.join(outdir,"summary.tsv"),sep="\t",header=True,index=False)


    if args.distance:
        dico_distance_sim = parse_dict_distance(distance)
        dico_distance_identity = parse_dict_distance(dist)
        if outdir is sys.stdout:
            print("Distance Identity : \n")
            pd.DataFrame.from_dict(dico_distance_identity,orient="index").to_csv(outdir,sep="\t",header=True,index=True)
            print("\nDistance Similarity : \n")
            pd.DataFrame.from_dict(dico_distance_sim,orient="index").to_csv(outdir,sep="\t",header=True,index=True)
        else:        
            pd.DataFrame.from_dict(dico_distance_sim,orient="index").to_csv(os.path.join(outdir,"distance_sim.tsv"),sep="\t",header=True,index=True)
            pd.DataFrame.from_dict(dico_distance_identity,orient="index").to_csv(os.path.join(outdir,"distance_idt.tsv"),sep="\t",header=True,index=True)
        

    #     with open(os.path.join())
    #     for i in UPGMA(li_sim, ide):
    # #     #dico_dis_simi = parse_dict(distance)
    # # pd.DataFrame.from_dict(,orient="index").to_csv(os.path.join(outdir,"simi_distance.tsv"),sep="\t",header=True,index=True)
    # # #pd.DataFrame.from_dict(UPGMA(li_iden, ide),orient="index").to_csv(os.path.join(outdir,"ident_distance.tsv"),sep="\t",header=True,index=True)
    