#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Jan 2022
# Modified: May 2022
# Description: convert an IQ-TREE Nexus concordance factor tree to Newick
##########################


import sys
import argparse
from Bio import Phylo


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to convert an IQ-TREE Nexus concordance factor tree to Newick')


# add arguments to parse
parser.add_argument('-t', type = str, dest = 'tree_file', help = 'The IQ-TREE Nexus file to convert')
parser.add_argument('-o', type = str, dest = 'out_pre', help = 'The output prefix [default \"output\"]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

tree_file = args.tree_file
out_pre = args.out_pre

if not tree_file:
	parser.print_help(sys.stderr)
	sys.exit(1)

if not out_pre:
	out_pre = 'output'


# Load the tree
tree = Phylo.read(tree_file, 'nexus')


# process the treefile to relabel internal nodes to the site concordance factors
# a typical comment looks like:
# [&label="100",sCF="43.44",sCF/sDF1/sDF2="43.44/42.83/13.73",sCF_N="391.94",sCF_N/sDF1_N/sDF2_N="391.94/386.88/123.89",sDF1="42.83",sDF1_N="386.88",sDF2="13.73",sDF2_N="123.89",sN="902.7089"]
# or if doing gene concordance factor:
# [&gCF="55.95",gCF/gDF1/gDF2/gDFP="55.95/19.05/20.24/4.76",gCF_N="47",gCF_N/gDF1_N/gDF2_N/gDFP_N="47/16/17/4",gDF1="19.05",
# gDF1_N="16",gDF2="20.24",gDF2_N="17",gDFP="4.76",gDFP_N="4",gN="84",label="100"]
#   gCF: Gene concordance factor (=gCF_N/gN %)
#   gCF_N: Number of trees concordant with the branch
#   gDF1: Gene discordance factor for NNI-1 branch (=gDF1_N/gN %)
#   gDF1_N: Number of trees concordant with NNI-1 branch
#   gDF2: Gene discordance factor for NNI-2 branch (=gDF2_N/gN %)
#   gDF2_N: Number of trees concordant with NNI-2 branch
#   gDFP: Gene discordance factor due to polyphyly (=gDFP_N/gN %)
#   gDFP_N: Number of trees decisive but discordant due to polyphyly
#   gN: Number of trees decisive for the branch
#   sCF: Site concordance factor averaged over sampled quartets (=sCF_N/sN %)
#   sCF_N: sCF in absolute number of sites
#   sDF1: Site discordance factor for alternative quartet 1 (=sDF1_N/sN %)
#   sDF1_N: sDF1 in absolute number of sites
#   sDF2: Site discordance factor for alternative quartet 2 (=sDF2_N/sN %)
#   sDF2_N: sDF2 in absolute number of sites
#   sN: Number of informative sites averaged over sampled quartets
#   Label: Existing branch label


# we want the second one (sCF/sDF1/sDF2)
scf_present = False
for node in tree.get_nonterminals():
	if node.comment:
		elements = str(node.comment).split(',')
		for element in elements:
			el_parts = element.split('=')
			if el_parts[0] == 'sCF/sDF1/sDF2':
				sitecf = el_parts[1].strip('"')
				scf_present = True
				break
		if scf_present:
			node.name = sitecf

# write the scf output tree
if scf_present:
	Phylo.write(tree, open(out_pre + '_scf.tre', 'w'), 'newick')


# for gcf, we want (gCF/gDF1/gDF2/gDFP)
gcf_present = False
for node in tree.get_nonterminals():
	if node.comment:
		elements = str(node.comment).split(',')
		for element in elements:
			el_parts = element.split('=')
			if el_parts[0] == 'gCF/gDF1/gDF2/gDFP':
				genecf = el_parts[1].strip('"')
				gcf_present = True
				break
		if gcf_present:
			node.name = genecf

# write the gcf output tree
if gcf_present:
	Phylo.write(tree, open(out_pre + '_gcf.tre', 'w'), 'newick')
