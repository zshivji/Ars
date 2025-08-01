import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

from itertools import combinations
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
import sys
from cluster_pos import cluster_pos

# get parsed hmm results
dir = sys.argv[1]
hits = pd.read_feather(f'../results/{dir}/hits.feather')

# save "contig" as col
hits['contig'] = hits['Hit'].str.split('_').str[:-1].str.join('_') 

# multi-index to cluster by genome, contig
hits.set_index(['GenomeID', 'contig'], inplace = True)
hits.sort_index(inplace = True)
hits.sort_values(by=['GenomeID', 'contig', 'Hit', 'E-value'], inplace=True) # drop duplicate hits, but keep most sig e-value
hits.drop_duplicates(inplace = True, subset=['Hit'], keep='first')

''' ArsB '''
# ArsB w/o ArsA --> keep all ArsB (later maybe add arsR,C) 
genomes_ArsB = hits[hits.Gene == 'arsB']

# filter for genomes to keep
genomes_ArsB.to_feather(f'../results/{dir}/arsB-only.feather')
genomes_ArsB.to_csv(f'../results/{dir}/arsB-only.csv')

''' ArsA '''
''' different than nif in that a bunch of ArsA can appear in same cluster --> makes order of filtering different than nif'''
# ArsA + ArsB --> filter for genome, contig with at least 2 unique genes (ArsAB)
filtered_arsA = hits.groupby(level=['GenomeID', 'contig']).filter(lambda x: x['Gene'].nunique() >= 2)

# make sure these 2 unique genes are not the same hit (i.e. not the same gene in reference genome)
filtered_arsA2 = filtered_arsA.groupby(level=['GenomeID', 'contig']).filter(lambda x: x['Hit'].nunique() >= 2)

genomes_ArsA = pd.DataFrame(columns = filtered_arsA2.columns) #ArsAB
# iterate through each genome and contig
for genome in filtered_arsA2.index.get_level_values(0).unique(): # iterate through each genome
    for contig in filtered_arsA2.loc[genome].index.get_level_values(0).unique(): # iterate through each contig

        tmp = filtered_arsA2.loc[(genome, contig)]

        # only keep numbers that have clusters >= 2 (clusters must be within 15 genes)
        pos_clusters = cluster_pos(tmp.Hit.unique())

        # for each cluster, find the best combination of genes (min e-value)
        for cl in pos_clusters:
            pos = [contig + '_' + str(p) for p in cl]
            no_pos = len(pos)
            
            # need at least 2 genes to continue
            if no_pos < 2:
                continue

            # only keep hits that are in the cluster
            tmp2 = tmp[tmp.Hit.isin(pos)].reset_index()

            # check if the cluster has at least 2 unique genes
            if tmp2.Gene.nunique() >= 2:
                genomes_ArsA = pd.concat([genomes_ArsA, tmp2])
        
# filter for genomes to keep
genomes_ArsA.set_index(['GenomeID', 'contig'], inplace = True)
genomes_ArsA.to_feather(f'../results/{dir}/arsA.feather')
genomes_ArsA.to_csv(f'../results/{dir}/arsA.csv')
