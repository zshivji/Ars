import pandas as pd
import os
import glob
from cluster_pos import cluster_pos

# grab checked ars
ars = pd.read_csv('../results/ars_rescheck_nofilt_07212025.csv')

#multi-index to cluster by genome, contig
ars.reset_index(inplace = True)
ars.set_index(['GenomeID', 'contig'], inplace = True)
ars.sort_index(inplace = True)
ars.drop_duplicates(inplace = True)

# iterate through each genome and contig
for genome in ars.index.get_level_values(0).unique(): # iterate through each genome
    for contig in ars.loc[genome].index.get_level_values(0).unique(): # iterate through each contig

        tmp = ars.loc[(genome, contig)]

        # get clusters
        pos_clusters = cluster_pos(tmp.Hit.unique(), 40)
        
        # get positions within +/-12 genes of center 
        for cl in pos_clusters:
            middle = cl[len(cl)//2] # center of cluster
            pos = [middle+num for num in range(-12,13) if middle+num > 0]
            acc = [contig + '_' + str(p) for p in pos] # get acc

            # write list items into a .txt file
            os.system("rm seqs.txt")
            with open("seqs.txt", "w") as f:
                for item in acc:
                    f.write(f"{item}\n")

            # save subsets as fasta
            try:
                file = glob.glob(f"../../GTDB/all_rep_proteins_aa/*/{genome[0]}_protein.faa")[0]
                os.system(f"seqtk subseq {file} seqs.txt > ../operon_org/{genome}_{contig}_operon.fasta")
            except:
                continue