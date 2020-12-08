#!/usr/bin/env python3

##############################################################################
### create_vector_merged.py
### Author: Tatiana Galochkina
### Takes as input the name for three encodings of the sequence (length N):
###  - pssm (*.aamtx_gaps: Nx21)
###  - aaindex (*.aaindex: Nx58)
###  - one-hot (*.onehot: Nx20)
### and the sliding window size. Concatenates the features 
### and adds the flag of terminus for each aa.
### Creates the output *.merged file of dimensions:
###  - N x (window_size (15) * total_number_of_fetatures (100))
################################################################################

import numpy as np
import argparse
import sys

def get_args():
    
    usage = "\nmerge_vector_features.py -i input_name -w window_size -o output_file"
    parser = argparse.ArgumentParser(usage = usage)
    parser.add_argument('-i', dest = "input_name", type = str, help = "Name\
    of the imput files")
    parser.add_argument('-w', dest = "window", type = int, default = 15, help = 
    "Sliding window size (15 by default)" )
    parser.add_argument('-o', dest = "output",
    type = str, help = "Output file to write merged vector" )
    parser.add_argument('-v', dest = "verbose", action = "store_true",
    help = "Verbose mode")
    args = parser.parse_args()
    
    return (args.input_name, args.window, args.output, args.verbose)


########################################
##-------------- MAIN-----------------##
########################################

# Reading input 
input_name, window, output, verbose = get_args()
aamtx = np.loadtxt(input_name+".aamtx_gaps")
aaind = np.loadtxt(input_name+".aaindex", skiprows=1)
onehot = np.loadtxt(input_name+".onehot", skiprows=1)

if (len(aamtx) != len(aaind)): 
    sys.stderr.write("aaindex and aamtx of different length!\n")
    sys.exit(1)

if (len(aamtx) != len(onehot)): 
    sys.stderr.write("onehot and aamtx of different length!\n")
    sys.exit(1)

if (len(onehot) != len(aaind)): 
    sys.stderr.write("aaindex and onehot of different length!\n")
    sys.exit(1)

# Total number of features
nfeatures = np.shape(aamtx)[1] + np.shape(aaind)[1] + np.shape(onehot)[1] + 1 # 99 in our case + 1 C-ter/N-ter flag
nres = len(aaind)
merged = np.zeros([nres, window*nfeatures])

# Indices for writing
start_aamtx = 0
start_aaind = np.shape(aamtx)[1]
start_onehot = np.shape(aamtx)[1] + np.shape(aaind)[1]

for i in range(nres):
    for j in range(-int(window/2), int(window/2)+1):
        if ((i+j) >= 0) and ((i+j) < nres): 
#            print(np.shape(merged[i, ((j+window-1)*nfeatures + start_aamtx):((j+window-1)*nfeatures + start_aamtx + np.shape(aamtx)[1])]))
#            print(((j + int(window/2))*nfeatures + start_aamtx))
#            print(((j + int(window/2))*nfeatures + start_aamtx + np.shape(aamtx)[1]))

            merged[i, ((j + int(window/2))*nfeatures + start_aamtx):((j + int(window/2))*nfeatures + start_aamtx + np.shape(aamtx)[1])] = aamtx[i+j, :]
            merged[i, ((j + int(window/2))*nfeatures + start_aaind):((j + int(window/2))*nfeatures + start_aaind + np.shape(aaind)[1])] = aaind[i+j, :]
            merged[i, ((j + int(window/2))*nfeatures + start_onehot):((j + int(window/2))*nfeatures + start_onehot + np.shape(onehot)[1])] = onehot[i+j, :]
        else:
            merged[i, (j + int(window/2))*nfeatures + start_onehot + np.shape(onehot)[1]] = 1

np.savetxt(output, merged, fmt = '%1.4f')

