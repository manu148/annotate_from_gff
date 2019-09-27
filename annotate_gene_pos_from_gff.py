import sys
import os
import gzip

gene_file = sys.argv[1]
output_file = sys.argv[2]
gff_file = "interim_GRCh37.p13_top_level_2017-01-13.gff3.gz"
chr_accession_file = "chr_accessions_GRCh37.p5"

gene_list = []
with open(gene_file) as gf:
	for line in gf:
		line = line.strip()
		gene_list.append(line)


accession_chr_map = {}
with open(chr_accession_file) as caf:
	for line in caf:
		line = line.strip()
		parts = line.split('\t')
		accession_chr_map[parts[1]] = parts[0]

gene_pos_map = {}
with gzip.open(gff_file) as gf:
	for line in gf:
		line = line.strip()
		parts = line.split('\t')
		if len(parts) > 2 and parts[2] == "gene" and parts[0] in accession_chr_map:
			annotations = parts[8].split(';')
			chrpos = accession_chr_map[parts[0]] + '\t' + parts[3] + '\t' + parts[4]
			for anno in annotations:
				if anno.startswith('gene='):
					gene_pos_map[anno.replace("gene=", '')] = chrpos
				if anno.startswith("gene_synonym="):
					synonyms = anno.replace("gene_synonym=", '')
					synonyms = synonyms.split(',')
					for synonym in synonyms:
						gene_pos_map[synonym] = chrpos

with open(output_file, 'w') as of:
	for gene in gene_list:
		if gene in gene_pos_map:
			of.write(gene_pos_map[gene] + '\t' + gene + '\n')
		else:
			print("Cannot find " + gene + " in gff")

