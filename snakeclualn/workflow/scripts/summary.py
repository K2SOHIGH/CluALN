import os
import sys
import pandas as pd

def clusummary(df):
    
    nof_clusters = len(list(df[df[0].duplicated()][0].unique()))
    nof_singletons = len(set(df[0])) - len(list(df[df[0].duplicated()][0].unique()))
    total = len(list(df[0]))
    t = df.groupby(0).count()
    return {
        "cluster": nof_clusters,
        "cluster (>10)": t[t[1]>=10].shape[0], 
        "input_in_cluster_>_10 (%)": t[t[1]>=10][1].sum() / total * 100, 
        "input_in_cluster (%)": t[t[1]>1][1].sum() / total * 100, 
        "singletons":nof_singletons,
        "singletons (%)":nof_singletons / total * 100,
        "total_seq":total,
    }

def per_clu_stats(clusterdir):
    empty = 0
    noempty = 0
    data = {}
    for cluster in os.listdir(clusterdir):
        table = os.path.join(clusterdir,cluster,"tables","summary.tsv")
        if (os.path.exists(table) and os.path.getsize(table) > 0):
            df = pd.read_csv(table,sep="\t",header=0,index_col=None)
            data[int(cluster.replace("cluster_",""))] = {
                "cluster_id":cluster,
                "cluster_rep":df.loc[0,"id1"],
                "#seq": len(set(list(df.id1) + list(df.id2))),                
                "align_len": df.len_aln.mean(),
                "mean_pid" : df.identity_percent.mean(),
                "mean_coverage" : df.coverage_1_2.mean(),
                "mean_reciprocal_coverage" : df.coverage_2_1.mean(),            
            }
            noempty += 1
        else:
            empty += 1
    return pd.DataFrame(data).T.reset_index().sort_values("index").set_index("index")

def per_clu_summary(df):
    return {
        "min_cluster_size":df["#seq"].min(),
        "max_cluster_size":df["#seq"].max(),        
        "mean_cluster_size":df["#seq"].mean(),
        "min_alignment_length":df["align_len"].min(),
        "max_alignment_length":df["align_len"].max(),        
        "mean_alignment_length":df["align_len"].mean(), 
        "min_mean_pid":df["mean_pid"].min(),
        "max_mean_pid":df["mean_pid"].max(),        
        "mean_mean_pid":df["mean_pid"].mean(), 
    }


resdir = str(snakemake.params.resdir)
clutable = os.path.join( resdir , "clustering/tables/clusters.tsv")
df = pd.read_csv(clutable , header=None ,sep="\t")
datas = clusummary(df)
df = pd.read_csv(clutable,sep="\t",header=None)

per_cluster_dir = os.path.join( resdir , "per_cluster")
per_clu_stats_df = None
if os.path.isdir(per_cluster_dir):
    per_clu_stats_df = per_clu_stats(per_cluster_dir)
    datas.update(per_clu_summary(per_clu_stats_df))
datas


outfile = str(snakemake.output)
with open(outfile,'w') as streamout:
    for i,j in datas.items():
        streamout.write("{} : {}\n".format(i,j))

    if per_clu_stats_df is not None:
        streamout.write("###############################\n")
        for i,j in per_clu_stats_df.T.to_dict().items():        
            for k,v in j.items():
                streamout.write("@cluster_{} - {} : {}\n".format(i,k,v))


