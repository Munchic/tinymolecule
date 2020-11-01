import os
import sys

sys.path.append("..")  # add top folder to path

import torch
import pandas as pd
import moses

import tinymolecule as tiny

# TO-DO add proper file path management

# load in download CCR4 data
data_dir = "/Users/Munchic/Developer/Capstone/tinymolecule/data/ccr4_ic50_meta.csv"
assays = pd.read_csv(data_dir)

# filter out assays with IC50
ic50_data = assays[assays["type"] == "IC50"]
ic50_nM_data = ic50_data[ic50_data["standard_units"] == "nM"]
ic50_nM_data.to_csv("/Users/Munchic/Developer/Capstone/tinymolecule/data/ccr4_ic50.csv")

# filter out valid SMILES codes
ccr4_ic50_smiles_df = pd.DataFrame({"SMILES": list(ic50_nM_data["canonical_smiles"])})
ccr4_ic50_smiles_df = ccr4_ic50_smiles_df[
    ccr4_ic50_smiles_df["SMILES"].apply(lambda x: isinstance(x, str))
]
ccr4_ic50_smiles_df.to_csv(
    "/Users/Munchic/Developer/Capstone/tinymolecule/data/ccr4_ic50_train.csv"
)

# create vocabulary
ccr4_ic50_smiles_vals = ccr4_ic50_smiles_df["SMILES"].values
vocab = moses.CharVocab.from_data(ccr4_ic50_smiles_vals)
torch.save(
    vocab, "/Users/Munchic/Developer/Capstone/tinymolecule/data/ccr4_ic50_vocab.pt"
)
