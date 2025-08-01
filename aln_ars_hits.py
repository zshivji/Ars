# align hits
from Bio import AlignIO
import pandas as pd
from Bio import SearchIO
from Bio import SeqIO
import os
import glob

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import numpy as np

# grab both archaea + bacteria hits
arsA_archaea = pd.read_feather('../results/archaea/arsA.feather')
# ars = ars_archaea.copy()
arsA_bacteria = pd.read_feather('../results/bacteria/arsA.feather')
ars = pd.concat([arsA_archaea, arsA_bacteria])

ars.reset_index(inplace = True)
ars.set_index(['GenomeID'], inplace = True)
ars['Seq'] = ''

# separate by annotation
arsA = ars[ars.Gene == 'arsA']
arsB = ars[ars.Gene == 'arsB']

gene_list = [arsA, arsB]
gene_names = ['arsA', 'arsB']

# get fasta sequences for each gene & export to fasta
for gene, name in zip(gene_list, gene_names):
    print(name)
    records = []
    for genome,hit in gene.iterrows():
        hit = hit.Hit

        file = glob.glob(f"../../GTDB/all_rep_proteins_aa/*/{genome}_protein.faa")[0]
        
        for result in SeqIO.parse(file, "fasta"):
            if result.id == hit:
                # store seq
                gene.loc[genome, 'Seq'] = str(result.seq)
                # convert to seqrecord
                record = SeqRecord(Seq(result.seq), id=genome, description=hit)
                records.append(record)
                # exit loop once sequence is found
                break
        
    # Write the records to a FASTA file
    with open("../results/" + name + ".fasta", "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

# align fasta files
for gene in ['arsA']:
    print(gene)
    num = eval(f"int({gene}.shape[0]/2000)+1") # how many splits
    os.system(f"seqtk split -n {num} ../results/fasta_splits/{gene}_split ../results/{gene}.fasta") # split fasta file
    for i in range(num):
        print(i+1)
        if i+1 < 10:
            os.system(f"seqtk subseq ../input-files/arsA_C113_172_422_all_seqs.fasta ref_seq.ids >> ../results/fasta_splits/{gene}_split.0000{i+1}.fa") # add reference sequences
            os.system(f"mafft --auto --quiet ../results/fasta_splits/{gene}_split.0000{i+1}.fa > ../results/fasta_splits/{gene}_split.0000{i+1}.aln")
        else:
            os.system(f"seqtk subseq ../input-files/arsA_C113_172_422_all_seqs.fasta ref_seq.ids >> ../results/fasta_splits/{gene}_split.000{i+1}.fa") # add reference sequences
            os.system(f"mafft --auto --quiet ../results/fasta_splits/{gene}_split.000{i+1}.fa > ../results/fasta_splits/{gene}_split.000{i+1}.aln")