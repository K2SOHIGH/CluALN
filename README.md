# SnakeClustering

Use mmseqs2 for sequence clustering. Then, each cluster is aligned using mafft and finaly MSAs are merged with mafft-merge utility. 


dependencies :
    - mmseqs
    - mafft
    - biopython
    
input [required] :
    - fasta file 
    - output directory path
    
input [optional] :
    - old_fasta : previous fasta file used for clustering
    - old_seqdb : -
    - old_cludb : -
    - coverage : - 
    - clumode : - 
    - verbose : [0..3]
output:
    -res_dir/clustering/tables/clusters.tsv 
    -res_dir/clustering/tables/identity.tsv 
    -res_dir/clustering/alignements/clusters.msa.aln
    
