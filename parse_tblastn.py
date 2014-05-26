#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import collections
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore -> std
# GPLIN_000706300 contig_3262     77.78   81      0       1       55      117     4938    4696    2e-33    126
# GPLIN_000706300 contig_3262     96.49   57      2       0       1       57      5731    5561    3e-19   85.9
# GPLIN_000707900 contig_98       100.00  64      0       0       15      78      529     720     5e-37    135
# GPLIN_000707900 contig_98       52.63   38      17      1       40      77      3815    3925    8e-05   41.6
# GPLIN_000707900 contig_98       64.00   25      9       0       15      39      3497    3571    3e-04   38.5
# GPLIN_000707900 contig_98       41.18   17      10      0       1       17      3363    3413    3e-04   22.3
# GPLIN_000707900 contig_98       100.00  14      0       0       1       14      329     370     0.019   33.9

def parse_blast_to_dict(filename):
	blast_dict = collections.defaultdict(dict) 
	with open(filename) as fh:
		for line in fh:
			temp = line.rstrip("\n").rsplit("\t")
			query, subject, sstart, send, strand = str(temp[0]), str(temp[1]), int(temp[8]), int(temp[9]), ''
			if sstart < send:
				strand = '+'
			else:
				strand = '-'   
			if query in blast_dict:
				if subject in blast_dict[query]:
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

def parse_assembly_to_dict(filename):
	assembly_dict = {}
	with open(filename) as fh:
		header = ''
		seq = ''
		for line in fh:
			line = line.rstrip("\n")
			if line.startswith(">"):
				header = line.lstrip(">")
				assembly_dict[header] = seq
        	else:
        		assembly_dict[header] += seq
        assembly_dict[header] = seq
	return assembly_dict

def print_blast_dict(blast_dict):
	for query in blast_dict:
		sys.stdout.write(query)
		for subject in blast_dict[query]:
			strand = blast_dict[query][subject][0]
			sstart = blast_dict[query][subject][1]
			send = blast_dict[query][subject][2]
			print ("\t" + subject + "\t" + strand + "\t" + str(sstart) + "\t" + str(send))


if __name__ == "__main__":
	blast_file = sys.argv[1]
	assembly_file = sys.argv[2]
	region = sys.argv[3]

	blast_dict = parse_blast_to_dict(blast_file)
	assembly_dict = parse_assembly_to_dict(assembly_file)
	print_blast_dict(blast_dict)
	# get_regions_from_assembly(assembly_dict, region)