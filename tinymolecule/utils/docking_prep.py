import os
import uuid
from pathlib import Path

import numpy as np
import pandas as pd

from openbabel import openbabel
import moses

PDBQT_DIR = Path("/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb")
SAMPLES_DIR = Path("/Users/Munchic/Developer/Capstone/moses/tinymolecule-out/samples")
GEN_DIR = Path("/Users/Munchic/Developer/Capstone/tinymolecule/data/gen")


def hash_smiles(smiles, hash_crop=8):
    namespace = uuid.NAMESPACE_URL
    hash_val = uuid.uuid5(uuid.NAMESPACE_URL, name=smiles)
    cropped_hash_val = str(hash_val)[:hash_crop]

    return cropped_hash_val


def filter_valid(
    gen_path=SAMPLES_DIR / "sample_100000_100e.csv",
    save_dir=GEN_DIR,
):
    gen = pd.read_csv(gen_path)

    valid = [
        (moses.utils.get_mol(molec) not in [None, np.nan]) for molec in gen["SMILES"]
    ]
    gen_valid = gen[valid].reset_index(drop=True)
    gen_valid["base64_id"] = gen_valid["SMILES"].apply(hash_smiles)

    save_dir = (
        GEN_DIR / "ccr5_ic50_train.csv"  # "valid_sample_1e5.csv"
    )  # TODO: automate to corresponding file name
    gen_valid.to_csv(save_dir)


def get_pdbqt(
    molec_path=GEN_DIR / "ccr5_ic50_train.csv", save_dir=PDBQT_DIR
):  # "valid_sample_1e5"
    valid_samples = pd.read_csv(molec_path)
    os.makedirs(save_dir / "train", exist_ok=True)  # "valid_sample_1e5"

    for i in valid_samples.index:
        print(
            f"converting molecule {valid_samples['base64_id'][i]} ({i}/{len(valid_samples)})"
        )

        # write temporary SMILES
        molec = valid_samples["SMILES"][i]
        base64_id = valid_samples["base64_id"][i]
        temp_file = save_dir / "train" / f"temp_{base64_id}.smi"  # "valid_sample_1e5"
        with open(temp_file, "w") as temp:
            temp.write(molec)

        # TODO: openbabel conversion settings, use FullConvert might be faster
        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats("smi", "pdbqt")
        conv.AddOption("gen3d", conv.GENOPTIONS)
        conv.AddOption("h", conv.GENOPTIONS)

        # convert temporary SMILES to PDBQT format
        pdbqt_file = save_dir / "train" / f"{base64_id}.pdbqt"  # "valid_sample_1e5"
        conv.OpenInAndOutFiles(str(temp_file), str(pdbqt_file))
        conv.Convert()
        os.remove(temp_file)


# filter_valid(
#     gen_path="/Users/Munchic/Developer/Capstone/tinymolecule/data/ccr5_ic50_train.csv"
# )
get_pdbqt(molec_path=GEN_DIR / "ccr5_ic50_train.csv")
