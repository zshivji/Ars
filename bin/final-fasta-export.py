import pandas as pd
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

ars = pd.read_csv('../results/ars_final_08042025.csv', index_col=[0,2])

# export fasta files
arsA = ars[(ars.Gene == 'arsA')]
arsB = ars[(ars.Gene == 'arsB')]
arsC = ars[(ars.Gene == 'arsC')]
arsR = ars[(ars.Gene == 'arsR')]
arsD = ars[(ars.Gene == 'arsD')]

gene_list = [arsA, arsB, arsC, arsR, arsD]
gene_names = ['arsA', 'arsB', 'arsC', 'arsR', 'arsD']

# get fasta sequences for each gene & export to fasta
for gene, name in zip(gene_list, gene_names):
    print(name)
    records = []

    for genome,contig in gene.iterrows():
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
    with open("../results/checked_" + name + "_08042025.fasta", "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")
