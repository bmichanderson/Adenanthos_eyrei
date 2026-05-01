#!/usr/bin/env python3

##########################
# Author: Benjamin Anderson
# Date: Oct 2022
# Description: substitute 10 consecutive Ns in fasta/multifasta files with another multiple (default 100)
# Note:	the script will overwrite the fasta file!
##########################


import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to substitute 10 consecutive Ns in ' +
	'fasta/multifasta files with another multiple (default 100).')


# add arguments to parse
parser.add_argument('fastas', type=str, help='The fasta files', nargs='*')
parser.add_argument('-m', type=int, dest='mult', help='The multiple of Ns to substitute')


# parse the command line
if len(sys.argv[1:]) == 0:              # if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()
input_fasta_files = args.fastas
multiple = args.mult


# read in the multiple, if present
if multiple:
	substitution = 'N' * multiple
else:
	substitution = 'N' * 100


# read in each fasta/multifasta and substitute, then overwrite
for fastafile in input_fasta_files:
	fasta_list = []
	with open(fastafile, 'r') as infile:
		fastas = SeqIO.parse(infile, 'fasta')
		for fasta in fastas:
			if 'NNNNNNNNNN' in fasta.seq:
				fasta.seq = Seq(str(fasta.seq).replace('NNNNNNNNNN', substitution))
			fasta_list.append(fasta)
	with open(fastafile, 'w') as outfile:
		for fasta in fasta_list:
			SeqIO.write(fasta, outfile, 'fasta')
