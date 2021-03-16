import os
import sys

sys.path.append("..")  # add top folder to path

import torch
import pandas as pd
import moses

# TODO add proper file path management

drugs_rm = {
    "Maraviroc": "Cc1nnc(C(C)C)n1[C@H]1C[C@@H]2CC[C@H](C1)N2CC[C@H](NC(=O)C1CCC(F)(F)CC1)c1ccccc1",
    "INCB-9471": "CCO[C@@H]1Cc2cc(C(F)(F)F)ccc2[C@H]1N1CCN(C2(C)CCN(C(=O)c3c(C)ncnc3C)CC2)C[C@@H]1C",
    "Aplaviroc": "CCCCN1C(=O)[C@@H]([C@H](O)C2CCCCC2)NC(=O)C12CCN(Cc1ccc(Oc3ccc(C(=O)O)cc3)cc1)CC2",
    "Vicriviroc": "COC[C@@H](c1ccc(C(F)(F)F)cc1)N1CCN(C2(C)CCN(C(=O)c3c(C)ncnc3C)CC2)C[C@@H]1C",
    "AZD5672": "CCN(C(=O)Cc1ccc(S(C)(=O)=O)cc1)C1CCN(CC[C@@H](c2cc(F)cc(F)c2)C2CCN(S(C)(=O)=O)CC2)CC1",
    "PF-04634817": "CO[C@@H]1COCC[C@@H]1N[C@@H]1CC[C@@](C(=O)N2C[C@@H]3C[C@H]2CN3c2cc(C(F)(F)F)ncn2)(C(C)C)C1",
    "Cenicriviroc": "CCCCOCCOc1ccc(-c2ccc3c(c2)/C=C(/C(=O)Nc2ccc([S@@+]([O-])Cc4cncn4CCC)cc2)CCCN3CC(C)C)cc1",
}

# load in download ccr5 data
data_dir = "/Users/Munchic/Developer/Capstone/tinymolecule/data/ccr5_ic50_meta.csv"
assays = pd.read_csv(data_dir)

# filter out assays with IC50
ic50_data = assays[assays["type"] == "IC50"]
ic50_nM_data = ic50_data[ic50_data["standard_units"] == "nM"]
ic50_nM_data.to_csv("/Users/Munchic/Developer/Capstone/tinymolecule/data/ccr5_ic50.csv")

# filter out valid SMILES codes
ccr5_ic50_smiles_df = pd.DataFrame({"SMILES": list(ic50_nM_data["canonical_smiles"])})
ccr5_ic50_smiles_df = ccr5_ic50_smiles_df[
    ccr5_ic50_smiles_df["SMILES"].apply(lambda x: isinstance(x, str))
]

# remove duplicates
ccr5_ic50_smiles_df.drop_duplicates(inplace=True)

# remove known drugs
ccr5_ic50_smiles_df = ccr5_ic50_smiles_df[
    ~ccr5_ic50_smiles_df["SMILES"].isin(drugs_rm.values())
]

ccr5_ic50_smiles_df.reset_index(inplace=True, drop=True)
ccr5_ic50_smiles_df.to_csv(
    "/Users/Munchic/Developer/Capstone/tinymolecule/data/ccr5_ic50_train.csv"
)

# create vocabulary
ccr5_ic50_smiles_vals = ccr5_ic50_smiles_df["SMILES"].values
vocab = moses.CharVocab.from_data(ccr5_ic50_smiles_vals)
torch.save(
    vocab, "/Users/Munchic/Developer/Capstone/tinymolecule/data/ccr5_ic50_vocab.pt"
)
