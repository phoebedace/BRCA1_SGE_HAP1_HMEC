#read_count_fractions.py

# python read_count_fractions.py [read_counts.txt] [read_counts2.txt]

import sys
import os
import subprocess
import operator

outfile_open = open(sys.argv[2],'w')
total_reads = 0
with open(sys.argv[1], 'r') as read_counts_file:
	for line in read_counts_file:
		data = line.strip().split(":")
		reads = int(data[1])
		sample = data[0]
		total_reads+= reads

with open(sys.argv[1], 'r') as read_counts_file:		
	for line in read_counts_file:
		data = line.strip().split(":")
		reads = int(data[1])
		sample = data[0]
		read_fraction = float(reads)/total_reads
		out_line = '\t'.join(data)+'\t'+str(read_fraction)
		outfile_open.write(out_line+'\n')
	outfile_open.write("total_reads"+'\t'+str(total_reads)+'\t'+"1.00")




