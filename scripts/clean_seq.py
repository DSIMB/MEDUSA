#! /usr/bin/env python
# -*- coding: ISO-8859-1 -*-

# AUTHOR GHOUZAM YASSINE
# MASTER 2 BIOINFORMATIQUE
# UNIVERSITE PARIS DIDEROT

# CE SCRIPT NETTOIE UNE SEQUENCE AU FORMAT FASTA, BRUT 

import sys
import os
import argparse

AA = ['A','C','E','D','G','F','I','H','K','M','L','N','Q','P','S','R','T','W','V','Y']


def header_and_seq(filename):
	""" Retourne le header et la sequence d'un fasta ou la sequence selon 
	le format (Fasta ou raw text)"""
	f = open(filename, 'r')
	lignes = f.readlines()
	f.close()
	if lignes[0][0] == '>':
		return (lignes[0],lignes[1:])
	else:
		return ('',lignes)


def clean_seq(seq, AA):
	""" Retourne la sequence en enlevant les caractéres non connus """
	clean_seq = []
	for i in xrange(len(seq)) :
		sub_seq = ''
		for j in xrange(len(seq[i])):
			if seq[i][j] in AA:
				sub_seq += seq[i][j]
		clean_seq.append(sub_seq)
	return clean_seq

def write_clean_seq(infile, outfile, AA):
	""" Ecrit une unique sequence dans un format Fasta """
	header , seq = header_and_seq(infile)
	clean_s = clean_seq(seq, AA)
	out = open(outfile , "w")
	out.write(header)
	for i in clean_s :
		out.write(i+"\n")
	out.close()


def get_args():
	""" Gestion des parametres d'entrée du programme
	"""
	
	usage = "clean_seq.py -in input_seq_file -out output_file"
	parser = argparse.ArgumentParser(usage=usage)
	parser.add_argument('-in', dest="infile", type = str,
	help = "Sequence file (input)")
	parser.add_argument('-out', dest="outfile", type = str,
        help = "Name of output file")

	args = parser.parse_args()
	
	return args.infile, args.outfile


def main():
	infile, outfile = get_args()
	write_clean_seq(infile, outfile, AA)

if __name__=="__main__":
    main()

