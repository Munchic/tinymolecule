import os
import uuid
import shlex
import subprocess

import random


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
