# ClusteringModule

Use mmseqs2 for sequence clustering. If --per-clu-msa is set, each cluster is aligned using mafft and MSAs produced are merged with mafft-merge utility. 


# DEPENDENCIES :
    - snakemake
    - pyyaml
    - pandas
    
# INPUT :
    - directory with one or more fasta file, fast afile or yaml file with path to fasta file(s).
    
# OUTPUT:
    -res_dir/<fasta_name>/clustering/tables/clusters.tsv 
    -res_dir/<fasta_name>/clustering/tables/identity.tsv 
    -res_dir/<fasta_name>/clustering/alignements/clusters.msa.aln

# INSTALLATION:
```bash
 git clone git@github.com:K2SOHIGH/ClusteringModule.git && cd ClusteringModule;
 mamba create -n clualn python=3 && conda activate clualn;
 pip3 install . ;
```

# USAGE : 
`clualn -i path/to/fasta -o path/to/resdir --per-clu-msa`
 
