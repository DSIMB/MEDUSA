#! /usr/bin/env python
# -*- coding: utf8 -*-


#===============================================================================
#===============================================================================

# NOTICE

'''

prend en argument un fichier .prof et renvoit les 10 classes d'accessibiltié
au solvant.


Pierre Edouard GUERIN                   7 avril 2014



'''



#===============================================================================
#===============================================================================

# MODULES

import os
import sys
import numpy
import argparse

#===============================================================================
#===============================================================================

#FONCTIONS

def acc_prof(famille):
	'''
	fonction qui prend en entrée le nom de la famille et ouvre le fichier .prof
	et retourne son acc et la reliability associée
	'''
	acc_res = []
	r_i = []
	fichier_entree = famille
	with open(fichier_entree,'r') as fichier_lecture:
		psa = fichier_lecture.readlines()
		for ligne in psa:
			if ligne[0] != '#' and ligne[0] != 'N':
				ligne = ligne.split("\t")
				#ajouter l'accessibilité au solvant pour chaque classe
				classe=12
				while classe < 22:							
					acc_res.append(float(ligne[classe])/100)
					classe+=1
				#ajouter la reliability pour chaque residu
				r_i.append((float(ligne[7])/10))					
	return acc_res, r_i


#===============================================================================
#===============================================================================

#Gestion des drapeaux et des arguments de la commande

parser = argparse.ArgumentParser()

parser.add_argument('-f', dest = "famille", type = str, help = "famille\
	HOMSTRAD prof acc prediction (file filtered_hssp.prof)")


#===============================================================================
#===============================================================================

#MAIN
args = parser.parse_args()
famille = args.famille

acc = acc_prof(famille)


ri = acc[1]
prof = acc[0]


i=0
j=-1
while i < len(prof):
	classe = 12
	while classe < 22:
		print "{0:.4f}".format(prof[i]),		
		classe+=1
		i+=1
	j+=1	
	print "{0:.4f}".format(ri[j]),
	#i+=1
	print "\n",
