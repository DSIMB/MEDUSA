#! /usr/bin/env python
# -*- coding: ISO-8859-1 -*-

# AUTHOR GHOUZAM YASSINE
# MASTER 2 BIOINFORMATIQUE
# UNIVERSITE PARIS DIDEROT

# CE PROGRAMME GENERE UN PROFIL DE FREQUENCES DE PB MODIFIEES
# (+ PSEUDO-COUNTS, WEIGHTING SCHEMES) A PARTIR D'UN ALIGNEMENT
# MULTIPLE AU FORMAT PIR OU AU FORMAT FASTA

import sys
try :
	import numpypy as np
except :
	import numpy as np

import argparse

AA = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o',
'p', '-']

BG_FREQ = {
'a' : 0.0384,
'b' : 0.0437,
'c' : 0.0827,
'd' : 0.1789,
'e' : 0.0221,
'f' : 0.0649,
'g' : 0.0129,
'h' : 0.0223,
'i' : 0.0155,
'j' : 0.0093,
'k' : 0.0530,
'l' : 0.0533,
'm' : 0.3005,
'n' : 0.0192,
'o' : 0.0254,
'p' : 0.0330,
}

Z_profile = [0.0419075, 0.0532328, 0.1107183, 0.2162724, 0.0283734, 0.07741762,
 0.01680114, 0.02728463, 0.01991704, 0.009494544, 0.05746235, 0.05579135, 
0.2035751, 0.01871807, 0.02587043, 0.03716331]


MTXCAR = "abcdefghijklmnop"
MTXIND = {}
for i in MTXCAR :
	MTXIND[i] = MTXCAR.index(i)


def ali_long_pir(filename):
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
	return cpt-2


def pir2list(filename):
	'''Parsage d'un alignment pir en liste de listes'''
	
	f = open(filename)
	l=f.readline()
	ali_longeur=ali_long_pir(filename)
	total=[]
	while l!="" :
		if l[0]=='>':
			l=f.readline()
			l=f.readline()
			un_ali = ""
			for i in xrange(ali_longeur):
				un_ali += l.replace("\n","")
				l=f.readline()
			total.append(un_ali[:-1])
	f.close()
	return total


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


def get_ali_format(filename):
	"""Trouve le format du fichier de 
	l'alignement (Deux formats supportés)"""
	msg = "ERROR : %s UNKNOWN OR UNSUPPORTED FORMAT FILE\
	\nSupported format : PIR and FASTA"%filename
	f = open(filename,'r')
	ligne1 = f.readline()
	ligne2 = f.readline()
	f.close()
	condition1 = ligne1[0] == '>' and ligne1[3] == ';'
	condition2 = ':' in ligne2
	if condition1 and condition2 :
		# Il s'agit d'un fichier PIR
		return "PIR-FORMAT"
	elif (((not condition1) and (not condition2))
                or (condition1 and (not condition2))):
		if '>' in ligne1:
			# Il s'agit d'un format fasta
			return "FASTA-FORMAT"
		else :
			sys.stderr.write(msg+'\n')
			sys.exit(1) 
	else :
		sys.stderr.write(msg+'\n')
		sys.exit(1)


def ali2list(filename):
	""" Stocke l'alignement du sous forme de liste """
	
	format_ali = get_ali_format(filename)
	
	if format_ali == "PIR-FORMAT":
		return pir2list(filename)
	return fasta2list(filename)


def occ_pos(ali,pos):
	"""renvoi l'occurence d'aa sur la position pos"""
	occ_aa = [0]*17 # liste des occ aa dans l'ordre
	for seq in ali:
		occ_aa[AA.index(seq[pos])] += 1
	return occ_aa


def freq_pos(ali,pos):
	"""renvoi les frequences d'aa sur la position pos"""
	occ = occ_pos(ali, pos)
	freq_aa = [0]*17
	for i in xrange(len(occ)):
		freq_aa[i] = occ[i]/float(len(ali))
	return freq_aa


def freq_prof(ali):
	nb_pos=len(ali[0])
	return np.array([freq_pos(ali,i) for i in xrange(nb_pos)])


def occ_prof(ali):
	nb_pos=len(ali[0])
	return np.array([occ_pos(ali,i) for i in xrange(nb_pos)])


def qia(a, i, mtx) :
	'''Retourne le score entre a et i dans mtx'''
	return mtx[MTXIND[a],MTXIND[i]]


def Rc1(mat_counts):
	""" Retourne le nb de résidus differents dans chaques colones""" 
	R = []
	for i in xrange(mat_counts.shape[0]):
		Rc = 0		
		for j in xrange(mat_counts.shape[1]-1):
			if mat_counts[i,j] != 0:
				Rc += 1
		R.append(Rc)
	return np.array(R)


def Rc2(mat_counts):
	""" Retourne le nb de résidus differents dans chaques colones""" 
	R = []
	for i in xrange(mat_counts.shape[0]):
		Rc = 0		
		for j in xrange(mat_counts.shape[1]):
			if mat_counts[i,j] != 0:
				Rc += 1
		R.append(Rc)
	return np.array(R)


def car_in_alipos(ali, pos, car):
	""" Retourne vrai si car est dans l'alignement ali à la position pos """
	for i in ali:
		if i[pos] == car:
			return True
	return False


def only_car_in_alipos(ali, pos, car):
	""" Retourne vrai si ali contient que seulement car à pos 
	de l'alignement"""
	for i in ali:
		if i[pos] != car:
			return False
	return True


def filter_ali_car(ali, car):
	""" Retourne les positions de l'alignement où car n'est pas compris """
	handled_pos = [i for i in xrange(len(ali[0])) 
	if not car_in_alipos(ali, i, car)]
	
	filtered_ali = []	
	for i in xrange(len(ali)):
		filtered_seq = ''		
		for j in handled_pos:
			filtered_seq += ali[i][j]
		filtered_ali.append(filtered_seq)
	return filtered_ali


def filter_ali_car_treshold(ali, car, treshold):
	""" Retire les positions de l'alignement où car 
	est présent à plus de treshold %"""
        handled_pos = [i for i in xrange(len(ali[0]))
        if freq_pos(ali, i)[AA.index(car)] < treshold ]
	
        filtered_ali = []
        for i in xrange(len(ali)):
                filtered_seq = ''
                for j in handled_pos:
                        filtered_seq += ali[i][j]
                filtered_ali.append(filtered_seq)
        return filtered_ali


def filter_ali_onlycar(ali, car):
	""" Retire les positions de l'alignement qui contiennent 
	exclusivement car et renvoie les positions retirées"""
	
	handled_pos = [i for i in xrange(len(ali[0])) 
	if not only_car_in_alipos(ali, i, car)]
	
	pos_removed = list(set(range(len(ali[0]))) - set(handled_pos))
	
	filtered_ali = []
	for i in xrange(len(ali)):
		filtered_seq = ''
		for j in handled_pos:
			filtered_seq += ali[i][j]
		filtered_ali.append(filtered_seq)
	return filtered_ali, handled_pos, pos_removed


def contains_non_AA(seq, AA):
	""" Retourne Vrai si la sequence contient un element ne faisant pas
	partie de AA"""
	for i in seq :
		if i not in AA:
			return True
	return False


def filter_ali_seq_gap(ali, seuil, AA):
	""" Retire de l'alignement les sequences qui ont plus de seuil % de 
	gaps ou qui ont des symboles ne faisant pas partie de la liste AA"""
	n = len(ali[0])
	filtered = []
	handled_seq = []
	for i in xrange(len(ali)) :
		if ( (ali[i].count("-")/(n*1.0)) < seuil and not contains_non_AA(ali[i],AA) and len(ali[i])==n):
			filtered.append(ali[i])
			handled_seq.append(i)
	
	removed_seq = list(set(range(len(ali))) - set(handled_seq))
	
	return filtered, handled_seq, removed_seq


def substitute_Z(ali):
	""" Substitue les caractéres Z par des gaps """
	return [i.replace('Z','-') for i in ali]


def bca(Bc, a, Nc, raw_occ_c, mtx_freq):
	A = 0	
	for i in MTXCAR:
		nci = raw_occ_c[AA.index(i)]
		qia_value = qia(a, i, mtx_freq)
		Qi = BG_FREQ[i]
		A += (nci/Nc) * (qia_value/Qi)
	return Bc * A


def occ_weighted_pos(ali, pos, Norm_weights):
	"""renvoi l'occurence d'aa sur la position pos ponderée par w"""
	occ_aa = [0]*17
	for i in xrange(len(ali)):
		occ_aa[AA.index(ali[i][pos])] += Norm_weights[i]
	return occ_aa


def occ_prof_weighted(ali, Norm_weights):
	nb_pos=len(ali[0])
	return (np.array([occ_weighted_pos(ali,i, Norm_weights)
	for i in xrange(nb_pos)]))


def read_profile(filename):
	"""Lit un profil à partir d'un fichier"""
	fichier = open(filename,"r")    # pointeur sur 1er ligne fichier
	return np.array([line[:-1].split() for line in fichier], dtype=float)


def get_args():
	
	usage = "\nPB2freq_pypy.py -al alignment_file -m matrix -gts gap_treshold_seq\
	(default 70 %%) -gtc gap_treshold_column (default 70 %%)"
	parser = argparse.ArgumentParser(usage = usage)
	parser.add_argument('-al', dest = "ali_file", type = str, help = "File\
	of multiple alignment (Pir of fasta format)")
	parser.add_argument('-m', dest = "matrix_file", type = str, help = 
	"Matrix of substitution for the pseudo-counts" )
	parser.add_argument('-gts', dest = "gap_treshold_seq", default = 90,
	type = float, help = "Sequences with a percentage of gaps superior than the\
	gap_treshold were not considered" )
	parser.add_argument('-gtc', dest = "gap_treshold_col", default = 70, 
	type = float, help = "Columns with a percentage of gaps superior than the\
        gap_treshold_col were not considered for the weighting scheme" )

	args = parser.parse_args()
	
	return (args.ali_file, args.matrix_file, args.gap_treshold_seq/100.0
	, args.gap_treshold_col/100.0)


######################################
#-------------- MAIN-----------------#
######################################

ali_file, matrix_file, gap_treshold_seq, gap_treshold_col = get_args()

# Lecture de l'alignment pour construire le profile 
alignts = ali2list(ali_file)

# On remplace les 'Z' par des gaps
alignts = substitute_Z(alignts)

if len(alignts) == 1 and list(set(alignts[0]))==["-"]:
	for i in alignts[0]:
		for j in Z_profile:
                        print "%.4f"%j,
                print ''
	sys.exit(0)

# On enléve les sequences qui ont plus de 70% de gaps et des caracteres inconnues
alignts, handled_seq, removed_seq = filter_ali_seq_gap(alignts, gap_treshold_seq, AA)

# On enléve les positions qui ne contiennent que des gaps
alignts, handled_pos, pos_removed = filter_ali_onlycar(alignts, '-')

# Lecture de la matrix de substitution pour les pseudo-counts 
mtxfreq = read_profile(matrix_file)

raw_occ = occ_prof(alignts)  # Occurences brutes
Rc = Rc1(raw_occ)  # Diversité à chaques colonnes 

alignts_g_f = filter_ali_car_treshold(alignts, '-', gap_treshold_col)  # Alignement filtré
filtered_occ = occ_prof(alignts_g_f) # occurences
Rc_without_g = Rc1(filtered_occ) # Diversité dans l'alignement filtré

#------------------------------------------------------#
# Weighting scheme de J.G Henikoff et S. Henikoff 1994 #
#------------------------------------------------------#
Norm_weights = []
nb_pos = len(alignts_g_f[0])

for i in alignts_g_f:
	Total = 0
	for j in xrange(len(i)):
		Total += 1.0/( Rc_without_g[j] * filtered_occ[j,AA.index(i[j])])
	Norm_weights.append(Total/nb_pos)

occ_wght = occ_prof_weighted(alignts, Norm_weights)  # Occurences pondérées

#---------------------------------------------------#
# Pseudo-counts de J.G Henikoff et S. Henikoff 1996 #
#---------------------------------------------------#
m = 5
Bcs = m * Rc  # Nombre totale de pseudo-counts à chaque colonnes

pseudo_occ = np.zeros((raw_occ.shape[0],16))

# Chaques colonnes de l'alignment
for i in xrange(occ_wght.shape[0]):
	
	# Nombre de counts totale dans la colonne i
	Nc = occ_wght[i,:-1].sum()/1.0
	# Nombre totale de pseudo-counts dans la colonne i
	Bc = Bcs[i]
	
	# Pour chaques residus
	for j in xrange(len(AA)-1): 
		nca = occ_wght[i,j]  # Nb de counts de l'aa j dans la colonne i
		# Calcul du pseudo-count pour chaque residus de la colonne
		bca_value = bca(Bc, AA[j], Nc, occ_wght[i,:-1], mtxfreq)
		pca =( ( ( Nc/(Nc + Bc) ) * (nca/Nc) )
		+( ( Bc/(Nc + Bc) ) * (bca_value/Bc) ) ) 
		pseudo_occ[i,j] = pca


### Affichage des frequences modifiées ###
# On check si des positions ont été retirées
if pos_removed != []:
	sys.stderr.write("Warning : %s profil generation %d "%(ali_file,len(pos_removed)) +
	"POSITIONS were removed \n")

if removed_seq != []:
	sys.stderr.write("Warning : %s profil generation %d "%(ali_file,len(removed_seq)) +
        "SEQUENCES were removed \n")

pos = 0
for i in xrange(len(handled_pos)+len(pos_removed)):
	# Les positions enlevées correspondent aux extrémitées du à la présence de 'Z'
	if i in pos_removed :
		for j in Z_profile:
			print "%.4f"%j,
		print ''
	
	else :
		for j in xrange(pseudo_occ.shape[1]):
			print "%.4f"%pseudo_occ[pos,j],
		print ''
		pos += 1

