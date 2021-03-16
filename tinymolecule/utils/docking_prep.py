import os
import uuid
from pathlib import Path

import numpy as np
import pandas as pd

from openbabel import openbabel
import moses

PDBQT_DIR = Path("/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb")
SAMPLES_DIR = Path("/Users/Munchic/Developer/Capstone/tinymolecule/data/samples")
GEN_DIR = Path("/Users/Munchic/Developer/Capstone/tinymolecule/data/gen")


def change_file_ext(filename, ext=None):
    ext_pos = filename.find(".")
    new_filename = filename[:ext_pos]
    if ext:
        new_filename += f".{ext}"

    return new_filename


def hash_smiles(smiles, hash_crop=8):
    namespace = uuid.NAMESPACE_URL
    hash_val = uuid.uuid5(uuid.NAMESPACE_URL, name=smiles)
    cropped_hash_val = str(hash_val)[:hash_crop]

    return cropped_hash_val


def filter_valid(
    samples_path=SAMPLES_DIR / "sample_100000_100e.csv",
    save_dir=GEN_DIR,
):
    gen = pd.read_csv(samples_path)

    valid = [
        (moses.utils.get_mol(molec) not in [None, np.nan]) for molec in gen["SMILES"]
    ]
    gen_valid = gen[valid].reset_index(drop=True)
    gen_valid["uuid"] = gen_valid["SMILES"].apply(hash_smiles)

    save_dir = (
        GEN_DIR / "ccr5_ic50_train.csv"  # "valid_sample_1e5.csv"
    )  # TODO: automate to corresponding file name
    gen_valid.to_csv(save_dir)


def get_pdbqt(
    molec_path=GEN_DIR / "ccr5_ic50_train.csv", save_dir=PDBQT_DIR
):  # "valid_sample_1e5"
    valid_samples = pd.read_csv(molec_path)
    save_folder = molec_path.stem
    os.makedirs(save_dir / save_folder, exist_ok=True)  # "valid_sample_1e5"

    for i in valid_samples.index:
        print(
            f"converting molecule {valid_samples['uuid'][i]} ({i + 1}/{len(valid_samples)})"
        )

        # write temporary SMILES
        molec = valid_samples["SMILES"][i]
        uuid = valid_samples["uuid"][i]
        temp_file = save_dir / save_folder / f"temp_{uuid}.smi"  # "valid_sample_1e5"
        with open(temp_file, "w") as temp:
            temp.write(molec)

        # TODO: openbabel conversion settings, use FullConvert might be faster
        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats("smi", "pdbqt")
        conv.AddOption("gen3d", conv.GENOPTIONS)
        conv.AddOption("h", conv.GENOPTIONS)

        # convert temporary SMILES to PDBQT format
        pdbqt_file = save_dir / save_folder / f"{uuid}.pdbqt"  # "valid_sample_1e5"
        conv.OpenInAndOutFiles(str(temp_file), str(pdbqt_file))
        conv.Convert()
        os.remove(temp_file)


# # get valid generated molecules
# filter_valid(samples_path=SAMPLES_DIR / "ccr5_valid_45e.csv")
# get_pdbqt(molec_path=GEN_DIR / "ccr5_valid_45e.csv")

# get valid training molecules
filter_valid(
    samples_path="/Users/Munchic/Developer/Capstone/tinymolecule/data/ccr5_train.csv",
)
get_pdbqt(molec_path=GEN_DIR / "ccr5_train.csv")