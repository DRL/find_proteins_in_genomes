#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import collections
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore -> std
# GPLIN_000706300 contig_3262     77.78   81      0       1       55      117     4938    4696    2e-33    126
# GPLIN_000706300 contig_3262     96.49   57      2       0       1       57      5731    5561    3e-19   85.9
# GPLIN_000707900 contig_98       100.00  64      0       0       15      78      529     720     5e-37    135
# GPLIN_000707900 contig_98       52.63   38      17      1       40      77      3815    3925    8e-05   41.6
# GPLIN_000707900 contig_98       64.00   25      9       0       15      39      3497    3571    3e-04   38.5
# GPLIN_000707900 contig_98       41.18   17      10      0       1       17      3363    3413    3e-04   22.3
# GPLIN_000707900 contig_98       100.00  14      0       0       1       14      329     370     0.019   33.9

def parse_blast_to_dict(blast_file):
	blast_dict = collections.defaultdict(dict) 
	with open(blast_file) as fh:
		for line in fh:
			temp = line.rstrip("\n").rsplit("\t")
			query, subject, sstart, send, strand = str(temp[0]), str(temp[1]), int(temp[8]), int(temp[9]), ''
			if sstart < send:
				strand = '+'
			else:
				strand = '-'   
			if query in blast_dict:
				if subject in blast_dict[query]:
					if strand == blast_dict[query][subject][0]:
						if strand == '+':
							if blast_dict[query][subject][2] < send:
								blast_dict[query][subject][2] = send
							if blast_dict[query][subject][1] > sstart:
								blast_dict[query][subject][1] = sstart
						else:
							if blast_dict[query][subject][2] > send:
								blast_dict[query][subject][2] = send
							if blast_dict[query][subject][1] < sstart:
								blast_dict[query][subject][1] = sstart
				else:
					blast_dict[query][subject] = [strand, sstart, send]
			else:  
				blast_dict[query][subject] = [strand, sstart, send]
	return blast_dict

def parse_fasta_to_dict(fasta_file):
	data_dict = {}
	with open(fasta_file) as fh:
		header = ''
		seq = ''
		for line in fh:
			line = line.rstrip("\n")
			if line.startswith(">"):
				if header == '':
					pass
				else:
					data_dict[header] = seq
				seq = ''
				header = line.split(' ')[0][1:]
			else:
				seq += line
		data_dict[header] = seq
	return data_dict

def print_blast_dict(blast_dict):
	for query in blast_dict:
		sys.stdout.write(query)
		for subject in blast_dict[query]:
			strand = blast_dict[query][subject][0]
			sstart = blast_dict[query][subject][1]
			send = blast_dict[query][subject][2]
			#print ("\t" + subject + "\t" + strand + "\t" + str(sstart) + "\t" + str(send))
			if strand == '+':
				print ("\t" + subject + "\t" + str(send - sstart))
			else:
				print ("\t" + subject + "\t" + str(sstart - send))

def print_regions(assembly_file, assembly_dict, blast_dict, region):
	for query in blast_dict:
		for subject in blast_dict[query]:
			sequence = assembly_dict[subject]
			strand = blast_dict[query][subject][0]
			#print query + ' ' + subject + ' ' + strand + ' hit ' + str(blast_dict[query][subject][1]) + ' ' + str(blast_dict[query][subject][2])
			#print sequence
			if (strand == '+'):
				start_region = blast_dict[query][subject][1] - region
				end_region = blast_dict[query][subject][2] + region
			else:
				start_region = blast_dict[query][subject][2] - region
				end_region = blast_dict[query][subject][1] + region
			if (start_region < 1):
				start_region = 0 
			if (end_region > len(sequence)):
				end_region = len(sequence)
			sequence_region = sequence[start_region:end_region]

			#print query + ' ' + subject + ' ' + str(start_region) + ' ' + str(end_region)
			#print sequence_region
			#filename = os.path.dirname(assembly_file) + '/' + assembly_file + '.' + query + '.' + subject + '.fa'
			filename = query + '.' + assembly_file + '.' + subject + '.' + strand + '.fa'
			fh = open(filename, "w")
			fh.write(">" + query + '.' + assembly_file + '.' + subject + '.' + str(region) + "\n" + sequence_region + "\n")

def print_proteins (protein_dict, blast_dict):
	print protein_dict
	for query in blast_dict:
		print query
		filename = query + '.' + assembly_file
		fh = open(filename, "w")
		fh.write(">" + query + '.' + assembly_file + "\n" + protein_dict[query] )

if __name__ == "__main__":
	blast_file = sys.argv[1]
	protein_file = sys.argv[2]
	assembly_file= sys.argv[3]
	region = int(sys.argv[4])
	blast_dict = parse_blast_to_dict(blast_file)
	protein_dict = parse_fasta_to_dict(protein_file)
	assembly_dict = parse_fasta_to_dict(assembly_file)

	print_regions(assembly_file, assembly_dict, blast_dict, region)
	print_proteins(protein_dict, blast_dict)
