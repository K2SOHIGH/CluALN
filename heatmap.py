import seaborn as sns
import pandas as pd
import sys
import os
from matplotlib import pyplot as plt




matrix = sys.argv[1]
outmap = sys.argv[2]
dirname = os.path.abspath(os.path.dirname(outmap))

os.makedirs(dirname,exist_ok=True)
df = pd.read_csv(matrix,sep="\t",header=0,index_col=0)

fig = sns.clustermap(df,cmap="summer")


# plt.plot(fig)
plt.savefig(outmap)
plt.show()