#! /bin/python3.7
#!/usr/bin/env python

import os
import sys
import pyfastx
import math
import itertools

if len(sys.argv)<3:
	print ('Usage: fasta file outputfile')
	sys.exit()
else:
	input1=sys.argv[1]
	output1=sys.argv[2]
outputfile = open(output1, 'w')
fa = pyfastx.Fasta(input1, uppercase=True)

chr21 = fa['chr21']
bin_dict = {}

bin_count = int(round(len(chr21) / 1000, 2))
bin_size = 1000
pairs = []
bins = []

print("starting process, the processing time may vary depending on the input file")
for x in range(bin_count):
	start_value = x * bin_size
	seq = chr21[start_value:start_value + bin_size]
	gc_content = seq.gc_content

	if not math.isnan(gc_content) and gc_content < 50:
		for ext_gcc, ext_start in bins:
			if abs(ext_start - start_value) < 1000000 and (abs(ext_gcc - gc_content) * 200 / (ext_gcc + gc_content)) < 2:
				bin_dict[ext_start] = True
				bin_dict[start_value] = True
				pairs.append((ext_start, start_value))
		bins.append((gc_content, start_value))

outputfile.write("no of pairs: " + str(len(pairs)) + "\n")
outputfile.write("no of bins that fits criteria: " + str(len(bin_dict)) + "\n")
outputfile.write("coordinates:\n")
for x, y in pairs:
	outputfile.write("x:" + str(x) + " y:" + str(y) + "\n")

print("output can be seen in output file")
# for i in range(len(gc_bins)):
# 	n = gc_bins[i]
#
# 	print(i)
#
# 	for j in range(len(gc_bins)):
# 		m = gc_bins[j]
#
# 		if i == j:
# 			continue
#
# 		# distance between bin must be less than 1mb. each bin is 1kb
# 		if abs(i - j) < 1000:
# 			if abs(m - n) * 200 / (m + n) < 2:
# 				pairs.append((i, j))


