#! /usr/bin/env python
# -*- coding: ISO-8859-1 -*-

# AUTHOR GHOUZAM YASSINE
# MASTER 2 BIOINFORMATIQUE
# UNIVERSITE PARIS DIDEROT

# COMPARE THE PSI-BLAST OUTPUT TO THE MTX OUTPUT AND ADJUST THIS ONE

import sys
import argparse
import re
import copy


def zeros(shape):
	retval = []
	for x in range(shape[0]):
		retval.append([])
		for y in range(shape[1]):
			retval[-1].append(0)
	return retval


def match_score(alpha, beta,match_award,mismatch_penalty,gap_penalty):
	if alpha == beta:
		return match_award
	elif alpha == '-' or beta == '-':
		return gap_penalty
	else:
		return mismatch_penalty


def finalize(align1, align2,match_award,mismatch_penalty,gap_penalty):
	align1 = align1[::-1]    #reverse sequence 1
	align2 = align2[::-1]    #reverse sequence 2

	i,j = 0,0

	#calcuate identity, score and aligned sequeces
	symbol = ''
	found = 0
	score = 0
	identity = 0
	for i in range(0,len(align1)):
		# if two AAs are the same, then output the letter
		if align1[i] == align2[i]:
			symbol = symbol + align1[i]
			identity = identity + 1
			score += match_score(align1[i], align2[i],match_award,mismatch_penalty,gap_penalty)

		# if they are not identical and none of them is gap
		elif align1[i] != align2[i] and align1[i] != '-' and align2[i] != '-':
			score += match_score(align1[i], align2[i],match_award,mismatch_penalty,gap_penalty)
			symbol += ' '
			found = 0

		#if one of them is a gap, output a space
		elif align1[i] == '-' or align2[i] == '-':
			symbol += ' '
			score += gap_penalty

	identity = float(identity) / len(align1) * 100
	
	return align1,align2,identity


def needle(seq1, seq2, match_award,mismatch_penalty,gap_penalty):
	m, n = len(seq1), len(seq2)  # length of two sequences

	# Generate DP table and traceback path pointer matrix
	score = zeros((m+1, n+1))	  # the DP table

	# Calculate DP table
	for i in range(0, m + 1):
		score[i][0] = gap_penalty * i
	for j in range(0, n + 1):
		score[0][j] = gap_penalty * j
	for i in range(1, m + 1):
		for j in range(1, n + 1):
			match = score[i - 1][j - 1] + match_score(seq1[i-1], seq2[j-1],match_award,mismatch_penalty,gap_penalty)
			delete = score[i - 1][j] + gap_penalty
			insert = score[i][j - 1] + gap_penalty
			score[i][j] = max(match, delete, insert)

	# Traceback and compute the alignment 
	align1, align2 = '', ''
	i,j = m,n # start from the bottom right cell
	while i > 0 and j > 0: # end toching the top or the left edge
		score_current = score[i][j]
		score_diagonal = score[i-1][j-1]
		score_up = score[i][j-1]
		score_left = score[i-1][j]

		if score_current == score_diagonal + match_score(seq1[i-1], seq2[j-1],match_award,mismatch_penalty,gap_penalty):
			align1 += seq1[i-1]
			align2 += seq2[j-1]
			i -= 1
			j -= 1
		elif score_current == score_left + gap_penalty:
			align1 += seq1[i-1]
			align2 += '-'
			i -= 1
		elif score_current == score_up + gap_penalty:
			align1 += '-'
			align2 += seq2[j-1]
			j -= 1

	# Finish tracing up to the top left cell
	while i > 0:
		align1 += seq1[i-1]
		align2 += '-'
		i -= 1
	while j > 0:
		align1 += '-'
		align2 += seq2[j-1]
		j -= 1

	return finalize(align1, align2, match_award,mismatch_penalty,gap_penalty)


def get_limit_query_seq(hh_out):
	filin = open(hh_out, 'r')
	ends = []
	for line in filin:
		if line[:11] == 'Q Consensus':
			line = line[11:-1].split(' ')
			while line[0] == '':
				del line[0]
			ends.append(int(line[0])-2)
			ends.append(int(line[-2])-1)
	filin.close()
	return min(ends), max(ends)


def write_adjusted_mtx(mtx_file, pos_begin, pos_end):
	""" Write the MTX adjusted to sequence """
	f_mtx = open(mtx_file)
	mtx = f_mtx.readlines()
	f_mtx.close()
	
	if int(mtx[0][:-1]) != (pos_end-pos_begin)+1 :
		sys.stderr.write("Warning %s truncated from pos 1 to "%(mtx_file)+
		"%i and from %i to %i\n"%(pos_begin+1, pos_end+1, int(mtx[0][:-1])))
		seq = mtx[1].replace("\n","")
		scores = mtx[14:]
		
		mtx_adjusted = copy.deepcopy(mtx[:14])
		mtx_adjusted[1] = seq[pos_begin : pos_end+1]+'\n'
		mtx_adjusted[0] = str(len(mtx_adjusted[1].replace("\n","")))+'\n'
		mtx_adjusted.extend(scores[pos_begin : pos_end+1])
		
		out = open("test.mtx","w")
		for i in mtx_adjusted:
			out.write(i)
		out.close()


def write_adjusted_mtx_v2(mtx_file, seqtokeep_file,match_award,mismatch_penalty,gap_penalty):
	""" Write the MTX adjusted to sequence """
	f_mtx = open(mtx_file)
	mtx = f_mtx.readlines()
	f_mtx.close()
	
	seqtokeep = "".join(open(seqtokeep_file).readlines()[1:]).replace("\n","").replace(" ","")
	
	seq_mtx = mtx[1].strip().replace("\n","")
	
	if seqtokeep != seq_mtx:
		seqtokeep_aligned , seq_mtx_aligned , identity = needle(seqtokeep,seq_mtx,match_award,mismatch_penalty,gap_penalty)
		sys.stderr.write("Warning %s truncated \nIdentity:%d\nSeq to keep : %s\nis not equal to\nSeq mtx     : %s\nAligned:\n%s\n%s\n"%(mtx_file,identity,seqtokeep,seq_mtx_aligned,seqtokeep_aligned , seq_mtx_aligned))
		
		scores = mtx[14:]
		#sys.stderr.write("Size:%d\n"%(len(scores)))
		scores_adjust = []
		for i in xrange(len(seqtokeep_aligned)):
			if seq_mtx_aligned[i] != '-' and seqtokeep_aligned[i] != '-':
				sys.stderr.write("i:%d Scores: %s \n"%(i,scores[i]))
				scores_adjust.append(scores[i])
		
		mtx_adjusted = copy.deepcopy(mtx[:14])
		mtx_adjusted[1] = seqtokeep+'\n'
		mtx_adjusted[0] = str(len(seqtokeep))+'\n'
		mtx_adjusted.extend(scores_adjust)

		out = open(mtx_file,"w")
		for i in mtx_adjusted:
			out.write(i)
		out.close()


def get_args():
	usage = ("adjust_mtx.py -in blast.out -in2 mtx_file -n nround")
	parser = argparse.ArgumentParser(usage=usage)
	parser.add_argument('-in', dest = "blast_out", type = str,
		help = "Blast output file")
	parser.add_argument('-in2', dest = "mtx_file", type = str,
		help = "mtx file")
	parser.add_argument('-in3', dest = "seq_tokeep", type = str,
		help = "seq to keep")	

	args = parser.parse_args()
	return (args.blast_out, args.mtx_file, args.seq_tokeep)


def main():
	blast_out, mtx_file, seqtokeep = get_args()
	#hhr       mtx
	pos_begin , pos_end = get_limit_query_seq(blast_out)
	#write_adjusted_mtx(mtx_file, pos_begin , pos_end)

	match_award      = 10
	mismatch_penalty = -1000
	gap_penalty      = -10 # both for opening and extanding
	
	write_adjusted_mtx_v2(mtx_file, seqtokeep, match_award,mismatch_penalty,gap_penalty)

if __name__=="__main__":
	main()
