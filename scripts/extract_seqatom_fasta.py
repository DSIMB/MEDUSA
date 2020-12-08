#! /usr/bin/env python
# -*- coding: ISO-8859-1 -*-

# AUTHOR GHOUZAM YASSINE
# DSIMB
# UNIVERSITE PARIS DIDEROT



import sys
import copy
import os
import re
import argparse
import PDB


def extract_seqatom(pdbfile,outfile):
	pdb_name = outfile.split("/")[-1].split(".")[0]
	
	try:
		pdb = PDB.PDB(pdbfile)
	except:
		PDB.clean(pdbfile)
		pdb = PDB.PDB(pdbfile)
	
	seqatom = ">%s\n"%pdb_name
	seqatom += pdb.aaseq()
	"""for i in pdb:
		if i[0][0:4]=="ATOM":
			seqatom += (PDB.AA1seq[PDB.resType(i[0][17:20].strip().replace("\n",""))])
	"""
	fich = open(outfile,"w")
	fich.write(seqatom)
	fich.close()


def get_args():
    usage = ("extract_seqatom_fasta.py -i pdbfile -o outfile")
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('-i', dest = "infile",  type = str, help = "Input pdbfile")
    parser.add_argument('-o', dest = "outfile", type = str, help = "Output fasta seq")
    args = parser.parse_args()
    return (args.infile, args.outfile)

def main():
    infile, outfile = get_args()
    extract_seqatom(infile,outfile)


if __name__=="__main__":
    main()

