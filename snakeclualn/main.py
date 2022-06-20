import argparse
import sys
import os
import yaml

_SNAKEFILE =  os.path.join(os.path.dirname(__file__), 'workflow/Snakefile')





def get_args():
    parser = argparse.ArgumentParser(
            prog='snakeclualn',
            description='Cluster sequences , perform MSA per cluster and merge MSA'
        )

    parser.add_argument(
        '-i',
        '--input',
        dest = "clu_input",
        required = True,
        help="a fasta file"
        )
    
    parser.add_argument(
        '-o',
        '--output-directory',
        dest = "res_dir",
        default = "resClualn",
        help = "output directory"
    )
    
    parser.add_argument(
        '-m',
        '--per-clu-msa',
        action = "store_true",
        help = "if set, each cluster will be aligned using MAFFT and then merged"
    )

    parser.add_argument(
        '-c',
        '--coverage',
        dest = "coverage",
        default=0.80,
        help="coverage threshold for clustering"
    )

    parser.add_argument(
        '--cm',
        dest="covmode",
        choices = [0,1,2,3],
        type = int,
        default=0,
        help="msmeqs covmode"
    )


    parser.add_argument(
        '--clumode',
        dest = "clumode",
        type = int,
        choices = [0,1,2],
        default = 1,
        help="mmseqs cluster mode"
    )

    parser.add_argument(
        '--pid',
        dest = "pid",
        type = float,
        default = 0,
        help="sequence identity threshold for clustering"
    )
    
    parser.add_argument(
        '--olddb',
        dest = "oldDB",
        type = str,
        default = None,
        help="path to file with for each line file : <new fasta id> \t <path to old fasta>  \t <path to old seqDB> \t <path to old cluDb>\n "
    )

    parser.add_argument(
        '-v',
        dest = "verbose",
        type = int,
        default = 1,
        help="mmseqs verbose"
    )

    parser.add_argument(
        '--log',
        dest = "log",
        type = str,
        default = None,
        help="logfile"
    )   

    parser.add_argument(
        '-e',
        '--extension',
        default=".fa",
        help="sequence file extension if input is a directory"
    )

    parser.add_argument('--snakargs', dest='snakargs', type=str, default="-j1",
            help='snakmake arguments')
    
    args = parser.parse_args()
    return args

def main():
    args = get_args()

    os.makedirs(args.res_dir , exist_ok= True)


    CONFIG = {}
    
    for arg in vars(args):
        if arg == "snakargs":
            SNAKARGS = getattr(args, arg)
        else:            
            #c = arg + "=" + getattr(args, arg)
            if arg == "clu_input":
                CONFIG[arg] = os.path.abspath(getattr(args,arg))
            else:
                CONFIG[arg] = getattr(args,arg) #.append(c)

    yaml.dump(CONFIG, open( args.res_dir + "/config.yaml" , 'w'))
        
    cmd = """
        snakemake --snakefile {snakefile} --use-conda --configfile {config} {snakargs}
    """.format( 
        snakefile = _SNAKEFILE , 
        config = args.res_dir + "/config.yaml", 
        snakargs = SNAKARGS 
    )
    
    print("running : " + cmd + " ... ")
    os.system(cmd)



if __name__ == "__main__":
    sys.exit(main())    