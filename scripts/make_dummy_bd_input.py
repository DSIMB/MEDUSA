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
	filout = open('dummy_' + f.replace('.mfasta', ''))
	seq = ''
	with open(input_f, 'r') as filin:
		for line in filin:
			if line[0] == '>':
				if seq != '':
					seq = seq.replace('-', '')
					i = 0
					while i < len(seq):
						filout.write(seq[i:i+80]+'\n')
						i += 80
					filout.write(seq[i-80:]+'\n')
				filout.write(line)
				
	filout.close()
