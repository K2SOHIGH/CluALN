# ClusteringModule

Use mmseqs2 for sequence clustering. Then, each cluster is aligned using mafft and finaly MSAs are merged with mafft-merge utility. 


# DEPENDENCIES :
    - mmseqs
    - mafft
    - biopython
    
# INPUT :
    - directory wuth one or more fasta file
    
# PARAMS :
    - \*oldDB : tsv file with : <new fasta id> \t <path to old fasta>  \t <path to old seqDB> \t <path to old cluDb>\n
    - coverage : - 
    - clumode : - 
    - verbose : [0..3]
# OUTPUT:
    -res_dir/<fasta_name>/clustering/tables/clusters.tsv 
    -res_dir/<fasta_name>/clustering/tables/identity.tsv 
    -res_dir/<fasta_name>/clustering/alignements/clusters.msa.aln
    

# USAGE :
    - edit config file or overwrite default value with --config key=value
    
`snakemake --use-conda -j4 -n`
    
    
\*OLDDB :
    old fasta file will be compared to the new one. If both share some sequences then update clustering will be performed else new sequences will be added to existing clustering.
