import pandas as pd
from Bio import SearchIO
from Bio import SeqIO
from Bio import Phylo
import os
import re
import glob
import regex

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from multiprocessing import Pool
from functools import partial
import numpy as np

ars = pd.read_csv('../results/ars_rescheck_nofilt_07212025.csv', index_col=[0,2])
arsB_only = pd.read_csv('../results/arsB-only.csv')
arsB_only.set_index(['GenomeID', 'Hit'], inplace = True)

# export fasta files
arsA = ars[(ars.Gene == 'arsA')]
arsB = ars[(ars.Gene == 'arsB')]

gene_list = [arsA, arsB, arsB_only]
gene_names = ['arsA', 'arsB', 'arsB_only']

# get fasta sequences for each gene & export to fasta
for gene, name in zip(gene_list, gene_names):
    print(name)
    records = []

    for genome,contig in gene.iterrows():
        print(genome)
        file = glob.glob(f"../../GTDB/all_rep_proteins_aa/*/{genome[0]}_protein.faa")[0]
        for result in SeqIO.parse(file, "fasta"):
            if result.id == genome[-1]:
                # store seq
                gene.loc[genome, 'Seq'] = str(result.seq)
                # convert to seqrecord
                record = SeqRecord(result.seq, id=genome[0], description=genome[-1])
                records.append(record)
                # exit loop once sequence is found
                break
    print(len(records))    
    # Write the records to a FASTA file
    with open("../results/checked_" + name + "_07212025.fasta", "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")
