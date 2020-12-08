#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Charlotte Périn _ 29/12/2016

import os
import sys


if __name__ == '__main__':

	# Error management
	if len(sys.argv) != 2:
		sys.exit("ERROR: one argument is needed (name of input file (mfasta format))")
	input_f = sys.argv[1]

	liste = []
	seq = ''
	maxlen = 0
	with open(input_f, 'r') as filin:
		for line in filin:
			if line[0] == '>':
				if seq != '':
					seq = seq.replace('-', '.')
					liste.append([name, seq])
					name = line.strip('\n')
					if (len(name)>12):
						name = name[0:12]
					name=name.replace('tr', '')
					name=name.replace('sp', '')
					name=name.replace('|', '')
					seq = ''
				else:
					name = 'query'
			else:
				seq += line[:-1]
	seq = seq.replace('-', '.')
	liste.append([name, seq])

	j = 0
	filout = open(input_f.replace('mfasta', 'saf'), 'w')
	filout.write('# SAF (Simple Alignment Format)\n#')
	while j < len(liste[-1][1]):
		for hit in liste:
			filout.write("\n{:<22s}".format(hit[0]))
			for i in xrange(j, j+50):
				if i == len(hit[1]):
					break
				if i % 10 == 0:
					filout.write(' ')
				filout.write(hit[1][i])
		filout.write('\n')
		j += 50
	filout.close()

