#!/usr/bin/python3.6
# Calculate average posterior per site for one individual
# Ex: getAvrLAI.py bLge-028 /media/labgenoma4/DATAPART1/fjtrindade/aHMM/max_vcf/run_minDP08/chrX/ ../oge_chr.txt
# genotropics workshop 2023

## Libraries and arguments
import sys
import os
import pandas as pd

my_args = sys.argv
if len(my_args) != 4:
	print("Usage: getAvrLAI.py sample_id path_posterior chr_file")
	exit()
else:
	sample = my_args[1]
	path_posterior = my_args[2]
	chr_file = my_args[3]

	input = path_posterior + sample + ".posterior"

	if os.path.isfile(input):
		pass
	else:
		print("Files dont exist!")
		exit()

	output = path_posterior + sample + ".posterior.parent1.csv"
	output2 = path_posterior + sample + ".posterior.avr.csv"

## Re-calculate the posterior considering 2 probs
ind_posterior = pd.read_csv(input, sep='\t', header=0)
ind_posterior['avr_parent1'] = list(ind_posterior["2,0"]+(ind_posterior["1,1"]/2))

## Seting color and pos for R plotting
ind_posterior['pos'] = list(range(0,len(ind_posterior['chrom'])))

chr_col = []
with open(chr_file, 'r') as fi:
	n=0
	for ln in fi:
		chr = []
		chr.append(ln.strip())

		if n==0:
			chr.append("darkseagreen")
			n=1
		else:
			chr.append("darkslategray")
			n=0

		chr_col.append(chr)

colors = []
for i in range(0,len(ind_posterior['chrom'])):
	for c in range(0,len(chr_col)):
		if ind_posterior['chrom'][i] == chr_col[c][0]:
			colors.append(chr_col[c][1])

ind_posterior['color'] = colors

## Summary average
parent1 = ind_posterior["avr_parent1"].mean()

summary = []
summary.append(str(ind_posterior["2,0"].mean()))
summary.append(str(ind_posterior["0,2"].mean()))
summary.append(str(ind_posterior["1,1"].mean()))
summary.append(str(parent1))
summary.append(str(1-parent1))

## Creating output files
ind_posterior.drop(['2,0','1,1','0,2'], axis=1, inplace=True)
ind_posterior.to_csv(output, index=False, sep='\t')

with open(output2, 'w') as fo:
	header = ["parent1_avr","parent2_avr","heteroz_avr","parent1","parent2"]
	fo.write(','.join(header))
	fo.write('\n')
	fo.write(','.join(summary))
