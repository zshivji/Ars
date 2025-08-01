from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd
import numpy as np
import os
import glob

def check_gene(gene, ref_seq, important_residues, passing_score, p=False):

    alignment = AlignIO.read(f"{gene}.aln", "fasta")

    for record in alignment:
        if ref_seq in record.description:
            aln = record.seq
            break

    # map important cols in aln to ref_seq
    residue_to_alignment = {}
    residue_idx = 0  # index in original (ungapped) sequence
    col2res = {}

    for aln_idx, char in enumerate(aln):
        if char != '-':
            residue_idx += 1
            if residue_idx in important_residues:
                residue_to_alignment[residue_idx] = aln_idx
                col2res[aln_idx] = aln[aln_idx]
                
            # early stop if we've found everything
            if len(residue_to_alignment) == len(important_residues):
                break

    # print the corresponding residues in the original sequence
    if p:
        for residue, aln_idx in residue_to_alignment.items():
            print(f"Residue {residue} corresponds to alignment index {aln_idx}: {aln[aln_idx]}")

    # store in dataframe        
    acc = [result.id for result in alignment]
    seqs = [list(str(result.seq)) for result in alignment]
    hits = [result.description.split(' ')[-1] for result in alignment]

    pssm = pd.DataFrame(seqs, index = acc)

    pssm = pssm.iloc[:, list(residue_to_alignment.values())]
    pssm['hit'] = hits
    pssm['contig'] = pssm['hit'].str.split('_').str[:-1].str.join('_')
    
    # check if cols contain correct residue for function
    def check_res(row):
        score = 0
        for col in residue_to_alignment.values():
            if row[col] == aln[col]:
                score += 2
                if aln[col] == 'C': # higher weight for correct C
                    score += 1
            elif aln[col] == 'C': # if ref seq is C, must also be C
                continue
            elif row[col] != '-': # greater penalty for gap than incorrect residue
                score += 1
        if score >= passing_score:
            return score
        else:
            return np.nan
        #return score

    pssm['score'] = pssm.apply(check_res, axis = 1)
    pssm.dropna(subset = ['score'], inplace = True)
    # rename cols
    pssm.rename(columns = col2res, inplace = True)

    return pssm

print("running", flush=True)

# arsA
print("checking arsA", flush=True)

important_residues_arsA = [113, 172, 422] # CCC in e. coli

ref_seq_arsA = 'tr|A0A2A5MC03|A0A2A5MC03_9ENTR'

passing_score = 9

# initialize dataframe
arsA_checked = check_gene('../results/fasta_splits/arsA_split.00001', ref_seq_arsA, important_residues_arsA, passing_score, p=True)

for file in glob.glob(f'../results/fasta_splits/arsA_split.000*.aln'):
    if file == '../results/fasta_splits/arsA_split.00001.aln':
        continue
    new = check_gene(file[:-4], ref_seq_arsA, important_residues_arsA, passing_score)
    arsA_checked = pd.concat([arsA_checked, new])

arsA_checked.drop_duplicates(subset = ['hit'], inplace = True)
arsA_checked.set_index(['hit'], append= True, inplace = True)

print(str(len(arsA_checked)) + " arsA seqs")

# append updated annotation (based on conserved residue matching) to ars files

# grab both archaea + bacteria hits
arsA_archaea = pd.read_feather('../results/archaea/arsA.feather')
arsA_bacteria = pd.read_feather('../results/bacteria/arsA.feather')
ars = pd.concat([arsA_archaea, arsA_bacteria])

ars.reset_index(inplace = True)
ars.set_index(['GenomeID', 'Hit'], inplace = True)
ars['Seq'] = ''
ars['residue_match'] = ''
ars['backup_match'] = ''

ars.sort_index(inplace = True)

# update residue match column in ars df
for gene in 'A':
    for genome, cols in eval(f"ars{gene}_checked.iterrows()"):
        ars.loc[(genome[0], genome[1]), 'residue_match'] = "ars" + gene

# filter to get hits that passed residue matching
ars.reset_index(level = 'Hit', inplace = True)
ars.set_index('contig', append=True, inplace = True)

arsA = ars[(ars.residue_match == 'arsA') & (ars.Gene == 'arsA') & (ars['Alignment Length'] > 400)]
arsB = ars.loc[arsA.index]
arsB = arsB[arsB.Gene == 'arsB'] # only keep arsB hits with true arsA

ars = pd.concat([arsA, arsB])
ars.sort_index(inplace = True)

ars.to_csv('../results/ars_rescheck_nofilt_07212025.csv', index = False)