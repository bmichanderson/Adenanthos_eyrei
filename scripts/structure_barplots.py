#!/usr/bin/env python3

##########################
# Author: B.M. Anderson
# Date: Nov 2021
# Modified: Mar 2022, Apr 2022, Dec 2023 (cleaned up; made default colours), Mar 2025 (better sample matching)
# Description: create barplots from Q matrices from Structure-like runs and output from CLUMPAK
##########################


import sys
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to create barplots for Structure-like Q matrices')


# add arguments to parse
parser.add_argument('-o', type = str, dest = 'out_pre', help = 'The prefix for the output pdf [default \"output\"]')
parser.add_argument('-c', type = str, dest = 'col_file', help = 'The colours file with one colour per line, e.g. #87cffa or green; ' +
	'this must be at least as long as the number of K to be plotted')
parser.add_argument('-p', type = str, dest = 'pops_file', help = 'The populations file in a tab-delimited form ' +
	'"sampleID    pop_num", one per line; samples should be in the same order as in the Q matrix')
parser.add_argument('-q', type = str, dest = 'Q_file', help = 'The Q matrix file containing the whitespace-delimited assignment proportions ' +
	'without headers or any line information; i.e. a table with K columns and as many rows as samples')
parser.add_argument('-s', type = str, dest = 'sorting', help = 'An optional file with the order of samples desired, one per line ' +
	'with the same IDs as in the pops_file')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
out_pre = args.out_pre
col_file = args.col_file
pops_file = args.pops_file
Q_file = args.Q_file
sorting = args.sorting

if any([not pops_file, not Q_file]):
	parser.print_help(sys.stderr)
	sys.exit(1)
if not out_pre:
	out_pre = 'output'


# load the colours file and capture
draw_colours = []
if col_file:
	with open(col_file, 'r') as infile:
		for line in infile:
			draw_colours.append(line.rstrip())
else:
	draw_colours.append('#539bff')
	draw_colours.append('#00ffb6')
	draw_colours.append('#008d15')
	draw_colours.append('#aa00df')
	draw_colours.append('#ffc800')
	draw_colours.append('#f1003e')
	draw_colours.append('#997336')
	draw_colours.append('#00ffef')
	draw_colours.append('#282db9')
	draw_colours.append('#ff78f5')
	draw_colours.append('#f4f622')


# load the pops file, keeping the order
sample_df = pd.read_csv(pops_file, sep = '\t', header = None, dtype = {0: 'string'})
sample_df.rename(columns = {0: 'Sample', 1: 'Pop'},	inplace = True)


# load the Q matrix
Q_df = pd.read_csv(Q_file, delim_whitespace = True, header = None)


# add the Q values to the sample dataframe
num_apops = len(Q_df.columns)
for K in range(num_apops):
	sample_df.loc[:, 'Ancestral_pop' + str(K + 1)] = Q_df[K]


# check that there are enough colours to complete the graph
if len(draw_colours) < num_apops:
	print('\nNot enough colours to complete the graph!\n')
	parser.print_help(sys.stderr)
	sys.exit(1)


# sort by pop, or if a specific sorting is provided, use that
if sorting:	
	sorted_samples = []
	with open(sorting, 'r') as infile:
		for line in infile:
			sorted_samples.append(line.rstrip())
	current_samples = list(sample_df.Sample)
	index = []
	for sample in sorted_samples:
		if sample in current_samples:
			index.append(current_samples.index(sample))
	plot_df = sample_df.loc[index]
else:
	plot_df = sample_df.sort_values(['Pop', 'Sample'])


# create the barplot
fig, ax = plt.subplots()
ax = plot_df.plot.bar('Sample', plot_df.columns[range(2, 2 + num_apops)], stacked = True,
	color = draw_colours[: num_apops], xlabel = '', ylabel = '',
	width = 1, edgecolor = 'black', figsize = (30, 5),
	legend = None)
for spine in ax.spines:
    ax.spines[spine].set_visible(False)
plt.tick_params(axis = 'y', left = False, labelleft = False)
plt.tick_params(axis = 'x', bottom = False)
plt.savefig(out_pre + '.svg', bbox_inches = 'tight', format = 'svg')
