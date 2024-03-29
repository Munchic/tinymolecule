import os
import shlex
import subprocess
from pathlib import Path
from random import sample
import itertools

import numpy as np
import pandas as pd


LOGS_SUMMARY_FILE = "summary.csv"


def change_file_ext(filename, ext=None):
    ext_pos = filename.find(".")
    new_filename = filename[:ext_pos]
    if ext:
        new_filename += f".{ext}"

    return new_filename


def shell_command(command):
    FNULL = open(os.devnull, "w")
    subprocess.check_call(shlex.split(command), stdout=FNULL, stderr=subprocess.STDOUT)


def get_vina_header():
    colnames = list(
        map(
            lambda x: [
                f"affin_kcal_mol-1_{x}",
                f"best_dist_rmsd_lb_{x}",
                f"best_dist_rmsd_ub_{x}",
            ],
            range(1, 11),
        )
    )

    num_sorted = sorted(
        list(itertools.chain.from_iterable(colnames)),
        key=lambda x: int(x[-x[::-1].find("_") :]),
    )
    cols_sorted = np.array(sorted(num_sorted, key=lambda x: x[: -x[::-1].find("_")]))

    return cols_sorted


def create_empty_entry():
    filler = np.ones((10, 3)) * np.nan

    return filler


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


def dock(
    ligand_folder_path,
    config_path,
    pdb_out_path,
    subsample_perc=0.01,
    molecs_from_logs=None,
):
    """
    Parameters
    ----------
        molecs_from_logs: Path or None
            Specifies path to a log file from which to take corresponding molecules for docking

    """

    # prepare paths
    out_path = pdb_out_path  # "valid_sample_1e5_ccr2"  # "train"  # "valid_sample_1e5"  # TODO: automate path choice
    logs_path = out_path / "logs"
    os.makedirs(out_path, exist_ok=True)
    os.makedirs(logs_path, exist_ok=True)

    all_ligand_files = list(  # filter out molecules whose structure we already have
        set(os.listdir(ligand_folder_path)) - set(os.listdir(out_path))
    )
    if os.path.isfile(
        str(molecs_from_logs)
    ):  # if choose molecules from logs of another docking experiment
        ligs_to_dock = pd.read_csv(molecs_from_logs)["uuid"].values
        ligand_files = []
        for lig in all_ligand_files:
            if change_file_ext(lig) in ligs_to_dock:
                ligand_files.append(lig)
    elif subsample_perc != False or subsample_perc < 1:
        subsample_count = int(subsample_perc * len(all_ligand_files))
        print(f"subsampling {subsample_count} molecules to dock")
        ligand_files = sample(all_ligand_files, subsample_count)

    else:
        ligand_files = all_ligand_files

    # dock every molecule
    for i, ligand in enumerate(ligand_files):
        print(f"docking molecule {i + 1} out of {len(ligand_files)}")
        try:
            shell_command(
                f"vina --config {config_path} "
                # + f"--receptor {receptor_path} "
                + f"--ligand {ligand_folder_path / ligand} "
                + f"--out {out_path / ligand} "
                + f"--log {logs_path / change_file_ext(ligand, ext='txt')}"
            )
        except subprocess.CalledProcessError as e:
            print(f"subprocess error at molecule {ligand}, skipping...")

    generate_logs_table(logs_path)


def generate_logs_table(logs_path):
    all_logs = os.listdir(logs_path)
    logs_table = get_vina_header()  # get header of log entries
    log_ids = []

    for logfile in all_logs:
        if logfile == LOGS_SUMMARY_FILE:  # ignore summary file that we will create
            continue

        with open(logs_path / logfile) as _log:
            lines = _log.readlines()

        if len(lines) > 27:  # where the table starts
            all_positions = lines[
                -lines[::-1].index("\n") + 3 : -1
            ]  # binding positions log
        else:
            continue

        for i in range(len(all_positions)):
            all_positions[i] = list(map(float, all_positions[i].split()))

        log = np.array(all_positions)[:, 1:]
        if log.shape == ():
            log = create_empty_entry()
        elif log.shape[0] < 10:
            filler = np.ones((10 - log.shape[0], 3)) * np.nan
            log = np.vstack((log, filler))
        entry = log.flatten(order="F")

        logs_table = np.vstack((logs_table, entry))
        log_ids.append(change_file_ext(logfile))

    logs_df = pd.DataFrame(logs_table[1:], columns=logs_table[0], dtype=float)
    logs_df["uuid"] = log_ids
    logs_df.to_csv(logs_path / LOGS_SUMMARY_FILE)  # save

    return logs_df


# # CCR5 docking on generated
# dock(
#     ligand_folder_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb/ccr5_valid_45e"
#     ),
#     config_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/vina_config.txt"
#     ),
#     pdb_out_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr5_valid_45e"
#     ),
#     subsample_perc=1,
# )

# generate_logs_table(
#     Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr5_valid_45e/logs"
#     )
# )

# # CCR5 docking on train
# dock(
#     ligand_folder_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb/ccr5_train"
#     ),
#     config_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/vina_config.txt"
#     ),
#     pdb_out_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr5_train"
#     ),
#     subsample_perc=1,
# )

# generate_logs_table(
#     Path("/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr5_train/logs")
# )

# # CCR2 docking on generated
# dock(
#     ligand_folder_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb/ccr5_valid_45e"
#     ),
#     config_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/vina_config_ccr2.txt"
#     ),
#     pdb_out_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr2_valid_45e"
#     ),
#     subsample_perc=1,
# )

# generate_logs_table(
#     Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr2_valid_45e/logs"
#     )
# )

# # CCR5 Y3D docking on generated
# dock(
#     ligand_folder_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb/ccr5_valid_45e"
#     ),
#     config_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/vina_config_ccr5_y3d.txt"
#     ),
#     pdb_out_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr5_y3d_valid_45e"
#     ),
#     subsample_perc=1,
# )

# generate_logs_table(
#     Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr5_y3d_valid_45e/logs"
#     )
# )

# # CCR5 G301E docking on generated
# dock(
#     ligand_folder_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb/ccr5_valid_45e"
#     ),
#     config_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/vina_config_ccr5_g301e.txt"f
#     ),
#     pdb_out_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr5_g301e_valid_45e"
#     ),
#     subsample_perc=1,
# )

# generate_logs_table(
#     Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr5_g301e_valid_45e/logs"
#     )
# )

# # CCR1 docking on generated
# dock(
#     ligand_folder_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb/ccr5_valid_45e"
#     ),
#     config_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/vina_config_ccr1.txt"
#     ),
#     pdb_out_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr1_valid_45e"
#     ),
#     subsample_perc=1,
# )

# generate_logs_table(
#     Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr1_valid_45e/logs"
#     )
# )

# # CCR3 docking on generated
# dock(
#     ligand_folder_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb/ccr5_valid_45e"
#     ),
#     config_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/vina_config_ccr3.txt"
#     ),
#     pdb_out_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr3_valid_45e"
#     ),
#     subsample_perc=1,
# )

# generate_logs_table(
#     Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr3_valid_45e/logs"
#     )
# )

# ordering: CCR2, CCR5 Y3D, CCR5 G301E, CCR1, CCR3
# ordering #2: CCR5 G202E, CCR5 I253F, CCR9, CCR6

# # CCR5 G202E docking on generated
# dock(
#     ligand_folder_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb/ccr5_valid_45e"
#     ),
#     config_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/vina_config_ccr5_g202e.txt"
#     ),
#     pdb_out_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr5_g202e_valid_45e"
#     ),
#     subsample_perc=1,
# )

# generate_logs_table(
#     Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr5_g202e_valid_45e/logs"
#     )
# )

# # CCR5 I253F docking on generated
# dock(
#     ligand_folder_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb/ccr5_valid_45e"
#     ),
#     config_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/vina_config_ccr5_i253f.txt"
#     ),
#     pdb_out_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr5_i253f_valid_45e"
#     ),
#     subsample_perc=1,
# )

# generate_logs_table(
#     Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr5_i253f_valid_45e/logs"
#     )
# )

# # CCR9 docking on generated
# dock(
#     ligand_folder_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb/ccr5_valid_45e"
#     ),
#     config_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/vina_config_ccr9.txt"
#     ),
#     pdb_out_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr9_valid_45e"
#     ),
#     subsample_perc=1,
# )

# generate_logs_table(
#     Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr9_valid_45e/logs"
#     )
# )

# # CCR6 docking on generated
# dock(
#     ligand_folder_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb/ccr5_valid_45e"
#     ),
#     config_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/vina_config_ccr6.txt"
#     ),
#     pdb_out_path=Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr6_valid_45e"
#     ),
#     subsample_perc=1,
# )

# generate_logs_table(
#     Path(
#         "/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb_out/ccr6_valid_45e/logs"
#     )
# )