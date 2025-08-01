# plot operon-org of annotated genes
import requests
import pandas as pd
import glob
import re
from Bio import SeqIO
from pygenomeviz import GenomeViz
import sys

# organize microbeannotator results
arsAB = pd.read_csv('../results/ars_rescheck_nofilt_07212025.csv', index_col=[0,1])

gene_abv = {'K00532': 'hydA', 'K07486':'tns', 'K07485':'tns', 'K07487':'tns', 'K13819': 'arsU'} # store ko2gene
gene_data = pd.DataFrame(columns = ['genome', 'contig', 'query_id', 'gene', 'ko_number', 'start', 'end', 'orientation'])
#gene_data.set_index(['query_id'], inplace = True)

def ko2gene(ko):
    try:
        return gene_abv[ko]
    except:
        url = f"https://rest.kegg.jp/get/ko:{ko}"
        r = requests.get(url)
        for line in r.text.split("\n"):
            if line.startswith("SYMBOL"):
                gene = line.split()[-1]
                gene_abv[ko] = gene
                return gene

# grab annotation files
for file in glob.glob(f"../operon-org/microbeannotator/annotation_results/*.annot"):
    # convet ko_number to gene abv
    annot = pd.read_csv(file, sep = '\t', index_col=0)
    annot['gene'] = annot['ko_number'].apply(ko2gene)

    # grab genome and contig from file name
    s = r'([A-Z_]+[A-Z0-9.]+)_([A-Z_]*[A-Z0-9.]+)'
    genome = re.search(s, file).group(1)
    contig = re.search(s, file).group(2)

    # grab start, end, orientation of each gene
    input = glob.glob(f"../operon-org/input-fastas/{genome}_{contig}_operon.fasta")[0]
    for result in SeqIO.parse(input, "fasta"):
        start = int(result.description.split('# ')[1])
        end = int(result.description.split('# ')[2])
        orientation = int(result.description.split('# ')[3])
        if result.id in annot.index:
            annot.loc[result.id, 'start'] = start
            annot.loc[result.id, 'end'] = end
            annot.loc[result.id, 'orientation'] = orientation

    # add arsA, arsB hits to annot
    for row in arsAB.loc[(genome, contig)].iterrows():
        annot.loc[row[1].Hit, 'gene'] = row[1].Gene
        annot.loc[row[1].Hit, 'start'] = int(row[1].Location.split('-')[0])
        annot.loc[row[1].Hit, 'end'] = int(row[1].Location.split('-')[1])
        annot.loc[row[1].Hit, 'orientation'] = 1 # BUG: CANNOT KNOW ORIENTATION

    # set genome and contig
    annot['genome'] = genome
    annot['contig'] = contig

    annot.reset_index(inplace = True)
    gene_data = pd.concat([gene_data, annot[['genome', 'contig', 'query_id', 'gene', 'ko_number', 'start', 'end', 'orientation']]])
gene_data.to_csv('../operon-org/operon-org-plot-data.csv')

# plot
to_plot = sys.argv[0]

if to_plot:
    # open microbeannotator results
    gene_data = pd.read_csv('../operon-org/operon-org-plot-data.csv', index_col=[1,2])
    genome_list = gene_data.index.get_level_values(0).unique().to_list()

    gv = GenomeViz(feature_track_ratio=0.7)
    gv.set_scale_bar(labelsize=20)
    # gv.set_scale_xticks(labelsize=20)

    for genome in genome_list:
        genes = gene_data.loc[(genome)].gene.to_list()
        starts = gene_data.loc[(genome)].start.to_list()
        ends = gene_data.loc[(genome)].end.to_list()
        orientations = gene_data.loc[(genome)].orientation.to_list()

        track = gv.add_feature_track(genome, (int(min(starts)), int(max(ends))))
        for idx, gene in enumerate(genes):
            if gene == 'arsA':
                color = '#156082'
            elif gene == 'arsR':
                color = '#76C0AF'
            elif gene == 'arsC':
                color = '#98A2C2'
            elif gene == 'arsB':
                color = '#C85D2B'
            elif gene == 'arsD':
                color = 'orange'
            else:
                color = '#FFFFFF'
                
            track.add_feature(
                int(starts[idx]),
                int(ends[idx]),
                int(orientations[idx]),
                plotstyle='bigarrow',
                arrow_shaft_ratio=0.7,
                fc = color,
                lw = 1,
                label = gene,
                text_kws=dict(rotation=0, vpos="bottom", hpos="center"))
            
        #track.add_sublabel()#int(min(starts)), int(max(ends))), fontsize=20)

    gv.savefig("../operon-org/test.png")
