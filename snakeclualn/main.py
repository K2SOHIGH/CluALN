import argparse
import sys
import os
import logging
import yaml

from snakeclualn.workflow.scripts import log
logger = log.setlogger()

_SNAKEFILE =  os.path.join(os.path.dirname(__file__), 'workflow/Snakefile')


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        #logging.FileHandler("clustering.log"),
        logging.StreamHandler()
    ]
)

def get_args():
    parser = argparse.ArgumentParser(
            prog='clualn',
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
        '--per-clu-msa',
        action = "store_true",
        help = "if set, each cluster will be aligned using MAFFT"
    )
    
    parser.add_argument(
        '-m',
        '--merge',
        action = "store_true",
        help = "if set, cluster MSA will be merged [warning ... might be slow]"
    )

    parser.add_argument(
        '-c',
        '--coverage',
        dest = "coverage",
        default=0.80,
        help="coverage threshold for clustering [0:1] (default: 0.8)"
    )

    parser.add_argument(
        '--cm',
        dest="covmode",
        choices = [0,1,2,3],
        type = int,
        default=0,
        help="mmsmeqs coverage mode (default: 0)"
    )


    parser.add_argument(
        '--clumode',
        dest = "clumode",
        type = int,
        choices = [0,1,2],
        default = 0,
        help="mmseqs cluster mode (default: 0)"
    )

    parser.add_argument(
        '--pid',
        dest = "pid",
        type = float,
        default = 0,
        help="sequence identity threshold for clustering [0:1] (default: 0)"
    )

    parser.add_argument(
        '-v',
        dest = "verbose",
        type = int,
        default = 1,
        choices = [0,1,2,3],
        help="mmseqs verbose [0 -> nothing, 1 -> +errors, 2 -> +warnings, 3 -> +info] "
    )

    parser.add_argument(
        '--log',
        dest = "log",
        type = str,
        default = None,
        help="logfile"
    )   

    # parser.add_argument(
    #     '-e',
    #     '--extension',
    #     default=".fa",
    #     help="sequence file extension if input is a directory"
    # )

    parser.add_argument('--snakargs', dest='snakargs', type=str, default="-j10",
            help='snakmake arguments')
    
    args = parser.parse_args()
    return args

def main():


    args = get_args()



    cprefix = os.path.abspath(
        os.path.join(
            os.path.expanduser('~'),
            "snakemake/clustering-module"
        )
    )
    os.makedirs(cprefix,exist_ok=True)

    os.makedirs(args.res_dir , exist_ok= True)

    logging.info("snakemake will install conda environment in %s" % cprefix)


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
        snakemake --snakefile {snakefile} --rerun-triggers mtime --use-conda --configfile {config} --conda-prefix {cp} {snakargs}
    """.format( 
        snakefile = _SNAKEFILE ,
        cp = cprefix,  
        config = args.res_dir + "/config.yaml", 
        snakargs = SNAKARGS 
    )
    
    logging.info("running : " + cmd + " ... ")
    os.system(cmd)
    logging.info("clustering done.")
    logging.info('exit')


if __name__ == "__main__":
    sys.exit(main())    
