#!/usr/bin/env python3
"""
DOCKSTRINGS FFS
"""


import os
import sys
import math
import datetime
import argparse
import h5py
import pandas as pd
import numpy as np
import functools

from tensorflow.keras.models import model_from_json
from tensorflow.keras.utils import to_categorical
from tensorflow.keras.losses import CategoricalCrossentropy
from tensorflow.keras.optimizers import Adam

from itertools import product
from sklearn.metrics import matthews_corrcoef


# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="Path to the merge feature vector file.", type=str, required=True)
parser.add_argument("-o", "--output", help="Path to the output directory.", type=str, required=True)
parser.add_argument("-m", "--models", help="Path to the models directory.", type=str, required=True)
parser.add_argument("-f", "--fasta", help="Path to the fasta file.", type=str, required=True)
parser.add_argument("-S", "--strict", help="Strict prediction.", action="store_true", default=True)
parser.add_argument("-NS", "--nonstrict", help="Non strict prediction.", action="store_true", default=True)
parser.add_argument("-3", "--three", help="3 classes prediction.", action="store_true", default=True)
parser.add_argument("-5", "--five", help="5 classes prediction.", action="store_true", default=True)
parser.add_argument("-d", "--dim", metavar="D", type=int, nargs=2, help="Rows Columns", required=True)
args = parser.parse_args()

MODEL = args.models

#Check if the merge vector file is valid
merge = args.i
if not os.path.isfile(merge):
    sys.exit(f"{merge} does not exist.\n"
			  "Please enter a valid merge vector file.")

OUTPUT = args.output

fasta = args.fasta
if not os.path.isfile(fasta):
    sys.exit(f"{fasta} does not exist.\n"
			  "Please enter a valid fasta file.")

if not (args.strict or args.nonstrict or three or five):
    parser.error("No model selected.")

model_list = []
if args.strict:
    model_list.append("S")
if args.nonstrict:
    model_list.append("NS")
if args.three:
    model_list.append("3")
if args.five:
    model_list.append("5")


ROWS, COLS = args.dim


if __name__ == "__main__":

    #Get entry files path
    ENTRY = merge

    ########################### LOADING MODELS #################################
    if "S" in model_list:

        PATH_MODEL_JSON = os.path.join(MODEL, "S", "model.json")
        PATH_MODEL_H5 = os.path.join(MODEL, "S", "weights_S.h5")

        # Load JSON arch
        with open(PATH_MODEL_JSON, "r") as json_file:
            loaded_model_json = json_file.read()
        model_S = model_from_json(loaded_model_json)

        # Load weights from H5
        model_S.load_weights(PATH_MODEL_H5)

        # Compile model
        model_S.compile(loss=CategoricalCrossentropy(),
                        optimizer=Adam())

    if "NS" in model_list:

        PATH_MODEL_JSON = os.path.join(MODEL, "NS", "model.json")
        PATH_MODEL_H5 = os.path.join(MODEL, "NS", "weights_NS.h5")

        # Load JSON arch
        with open(PATH_MODEL_JSON, "r") as json_file:
            loaded_model_json = json_file.read()
        model_NS = model_from_json(loaded_model_json)

        # Load weights from H5
        model_NS.load_weights(PATH_MODEL_H5)

        # Compile model
        model_NS.compile(loss=CategoricalCrossentropy(),
                        optimizer=Adam())

    if "3" in model_list:

        PATH_MODEL_JSON = os.path.join(MODEL, "3", "model.json")
        PATH_MODEL_H5 = os.path.join(MODEL, "3", "weights_3.h5")

        # Load JSON arch
        with open(PATH_MODEL_JSON, "r") as json_file:
            loaded_model_json = json_file.read()
        model_3 = model_from_json(loaded_model_json)

        # Load weights from H5
        model_3.load_weights(PATH_MODEL_H5)

        # Compile model
        model_3.compile(loss=CategoricalCrossentropy(),
                        optimizer=Adam())

    if "5" in model_list:

        PATH_MODEL_JSON = os.path.join(MODEL, "5", "model.json")
        PATH_MODEL_H5 = os.path.join(MODEL, "5", "weights_5.h5")

        # Load JSON arch
        with open(PATH_MODEL_JSON, "r") as json_file:
            loaded_model_json = json_file.read()
        model_5 = model_from_json(loaded_model_json)

        # Load weights from H5
        model_5.load_weights(PATH_MODEL_H5)

        # Compile model
        model_5.compile(loss=CategoricalCrossentropy(),
                        optimizer=Adam())




    #Load X
    X = np.loadtxt(ENTRY)
    x_test = np.reshape(X, (len(X), ROWS, COLS))

    #Load fasta file
    fasta_file = np.loadtxt(fasta, dtype="str", skiprows=1)
    seq_AA = fasta_file.tolist()

    #Create output directory
    os.makedirs(OUTPUT, exist_ok=False)


    if "S" in model_list:
        y_pred = model_S.predict(x_test, verbose=0)
        y_pred_classes = np.argmax(y_pred, axis=1)
        prob_y_pred = np.amax(y_pred, axis=1)

        # Create a dataframe from zipped list
        zippedList =  list(zip(seq_AA, y_pred_classes, prob_y_pred, y_pred[:, 0], y_pred[:, 1]))
        df = pd.DataFrame(zippedList, columns = ["res", "S", "P_max", "P_0", "P_1"])
        df.to_csv(os.path.join(OUTPUT, "S_prediction.csv"), index=False, float_format="%.2f", sep="\t")


    if "NS" in model_list:
        y_pred = model_NS.predict(x_test, verbose=0)
        y_pred_classes = np.argmax(y_pred, axis=1)
        prob_y_pred = np.amax(y_pred, axis=1)

        # Create a dataframe from zipped list
        zippedList =  list(zip(seq_AA, y_pred_classes, prob_y_pred, y_pred[:, 0], y_pred[:, 1]))
        df = pd.DataFrame(zippedList, columns = ["res", "NS", "P_max", "P_0", "P_1"])
        df.to_csv(os.path.join(OUTPUT, "NS_prediction.csv"), index=False, float_format="%.2f", sep="\t")


    if "3" in model_list:
        y_pred = model_3.predict(x_test, verbose=0)
        y_pred_classes = np.argmax(y_pred, axis=1)
        prob_y_pred = np.amax(y_pred, axis=1)

        # Create a dataframe from zipped list
        zippedList =  list(zip(seq_AA, y_pred_classes, prob_y_pred, y_pred[:, 0], y_pred[:, 1], y_pred[:, 2]))
        df = pd.DataFrame(zippedList, columns = ["res", "pred_3", "P_max", "P_0", "P_1", "P_2"])
        df.to_csv(os.path.join(OUTPUT, "3_prediction.csv"), index=False, float_format="%.2f", sep="\t")


    if "5" in model_list:
        y_pred = model_5.predict(x_test, verbose=0)
        y_pred_classes = np.argmax(y_pred, axis=1)
        prob_y_pred = np.amax(y_pred, axis=1)

        # Create a dataframe from zipped list
        zippedList =  list(zip(seq_AA, y_pred_classes, prob_y_pred, y_pred[:, 0], y_pred[:, 1], y_pred[:, 2], y_pred[:, 3], y_pred[:, 4]))
        df = pd.DataFrame(zippedList, columns = ["res", "pred_5", "P_max", "P_0", "P_1", "P_2", "P_3", "P_4"])
        df.to_csv(os.path.join(OUTPUT, "5_prediction.csv"), index=False, float_format="%.2f", sep="\t")
