import os
import uuid
import shlex
import subprocess
import random

import pandas as pd


def change_file_ext(filename, ext=None):
    ext_pos = filename.find(".")
    new_filename = filename[:ext_pos]
    if ext:
        new_filename += f".{ext}"

    return new_filename


def hash_smiles(smiles, hash_crop=8, hash_seed=42):
    rd = random.Random()
    rd.seed(hash_seed)
    namespace = uuid.NAMESPACE_URL
    hash_val = uuid.uuid5(uuid.NAMESPACE_URL, name=smiles)
    cropped_hash_val = str(hash_val)[:hash_crop]

    return cropped_hash_val


def shell_command(command):
    FNULL = open(os.devnull, "w")
    subprocess.check_call(shlex.split(command), stdout=FNULL, stderr=subprocess.STDOUT)


def get_smiles_from_id(logs_path, samples_path):
    req_ids = pd.read_csv(logs_path)["uuid"].values  # requested IDs
    samples = pd.read_csv(samples_path)

    smi_corresp = samples[samples["uuid"].isin(req_ids)]  # corresponding SMILES codes
    smi_corresp.drop_duplicates(subset="SMILES", inplace=True)
    smi_corresp.reset_index(drop=True, inplace=True)

    return smi_corresp


def get_mean_baff_df(logs, prefix=""):
    all_baff = pd.concat([logs[f"affin_kcal_mol-1_{i}"] for i in range(1, 11)], axis=1)
    mean_baff = all_baff.mean(axis=1) * (-1)
    mean_baff.rename(f"{prefix}avg_affin_kcal_mol-1", inplace=True)
    mean_baff_df = pd.concat([logs["uuid"], mean_baff], axis=1)

    return mean_baff_df