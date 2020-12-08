#! /usr/bin/env python
# -*- coding: ISO-8859-1 -*-

# AUTHOR GHOUZAM YASSINE
# MASTER 2 BIOINFORMATIQUE
# UNIVERSITE PARIS DIDEROT

# COMPARE THE QUERY SEQUENCE TO THE PSI-BLAST OUTPUT AND ADJUST THE QUERY SEQUENCE

import sys
import argparse


def get_last_round(blast_out_file, n_round):
	""" Take in parameter one blast output, and get hits from the last 
	round """
	f = open(blast_out_file,"r")
	lines = f.readlines()
	f.close()
	
	if n_round==1:
		return lines
	
	id_search = []
	for r in xrange(n_round):
		for i in xrange(len(lines)):
			if lines[i][:-1]=="Results from round "+str(r+1) :
				id_search.append(i)

	return lines[max(id_search):]


def ali_long_fasta(filename):
        '''Longueur de l'alignement en lignes'''
        f = open(filename)
        l = f.readline()

        flag=False
        cpt=0
        while l!="":
                if flag and l[0]=='>':
                        f.close()
                        break
                if l[0]=='>':
                        flag=True
                cpt+=1
                l = f.readline()
        return cpt-1


def fasta2list(filename):
        '''Parsage d'un alignment multiple fasta en liste de listes'''

        f = open(filename)
        l = f.readline()
        ali_longeur = ali_long_fasta(filename)
        total = []
        while l != "" :
                if l[0] == '>':
                        l = f.readline()
                        un_ali = ""
                        for i in xrange(ali_longeur):
                                un_ali += l.replace("\n","")
                                l = f.readline()
                        total.append(un_ali)
        f.close()
        return total


def headers_alifasta(filename):
        """ Retourne les headers du fichier d'aligment pir"""
        f = open(filename)
        l = f.readline()
        ali_longeur = ali_long_fasta(filename)
        headers = []
        while l!="" :
                if l[0]=='>':
                        headers.append(l[:-1])
                        l=f.readline()
                        for i in xrange(ali_longeur):
                                l = f.readline()
	
        f.close()
        return headers


def split_by_n(seq, n):
        return [seq[i:i+n] for i in range(0, len(seq), n)]


def write_fasta(outname, seq, header):
        f = open(outname, "w")
        f.write("%s\n"%header)
        formated_seq = split_by_n(seq, 79)
	
        for i in formated_seq:
                f.write(i+'\n')
        f.close()


def write_fasta_query(mfasta_name, filename_query):
	
	ali = fasta2list(mfasta_name)
	query_seq = ali[0].replace("\n","").replace("-","")
	query_header = headers_alifasta(mfasta_name)[0].replace("\n","")
	write_fasta(filename_query, query_seq, query_header)


def get_args():
        usage = ("adjust_seq.py -f1 ali.mfasta -f2 query_filename")
        parser = argparse.ArgumentParser(usage=usage)
        parser.add_argument('-in', dest = "infile", type = str,
        help = "Multiple alignment fasta file")
        parser.add_argument('-out', dest = "query_name", type = str,
        help = "query_filename")	

        args = parser.parse_args()

        return (args.infile, args.query_name)


def main():
	infile, outfile = get_args()
	write_fasta_query(infile, outfile)

if __name__=="__main__":
	main()
