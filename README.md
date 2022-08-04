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
 mamba install pyyaml snakemake pandas
 pip3 install . ;
```

# USAGE : 
```bash
usage: clualn [-h] -i CLU_INPUT [-o RES_DIR] [-m] [-c COVERAGE] [--cm {0,1,2,3}] [--clumode {0,1,2}] [--pid PID] [-v VERBOSE] [--log LOG] [-e EXTENSION] [--snakargs SNAKARGS]

Cluster sequences , perform MSA per cluster and merge MSA

optional arguments:
  -h, --help            show this help message and exit
  -i CLU_INPUT, --input CLU_INPUT
                        a fasta file
  -o RES_DIR, --output-directory RES_DIR
                        output directory
  -m, --per-clu-msa     if set, each cluster will be aligned using MAFFT and then merged
  -c COVERAGE, --coverage COVERAGE
                        coverage threshold for clustering
  --cm {0,1,2,3}        msmeqs covmode
  --clumode {0,1,2}     mmseqs cluster mode
  --pid PID             sequence identity threshold for clustering
  -v VERBOSE            mmseqs verbose
  --log LOG             logfile
  -e EXTENSION, --extension EXTENSION
                        sequence file extension if input is a directory
  --snakargs SNAKARGS   snakmake arguments
```
 
