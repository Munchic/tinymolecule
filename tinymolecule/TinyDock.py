import os
from copy import deepcopy
from pathlib import Path
from typing import Optional, Union, List, Tuple
from random import sample

import yaml
import numpy as np
import pandas as pdf

from openbabel import openbabel
from tinymolecule.utils.helper_fn import change_file_ext, hash_smiles
from tinymolecule.utils.molecules import molecule_is_valid


class TinyDock:
    """
    Core class of the tinymolecule pipeline to dock generated molecules.
    Stores docked molecules objects and their metrics. Implements analyses
    showcased in the paper.
    """

    _DEFAULT_CONFIG = Path(__file__).parent / "config" / "default_config.yaml"

    def __init__(self, config=None):
        with open(self._DEFAULT_CONFIG) as config_file:
            self.config = (
                yaml.load(config_file, Loader=yaml.Loader) if not config else config
            )
        self.paths = self.config["data_paths"]
        self.data_path = Path(self.config["data_paths"]["data"])
        self.targets = self.config["targets"]

        self.ligands_train_dir = (
            self.data_path / self.paths["ligands_dir"] / self.paths["ligands_train_dir"]
        )
        self.ligands_gen_dir = (
            self.data_path / self.paths["ligands_dir"] / self.paths["ligands_gen_dir"]
        )
        self.which_ligands = "gen"

        self.tiny_params = self.config["tiny_params"]

    def prepare_molecules(self, smiles_csv_path: Union[string, Path] = None):
        """
        Creates pdbqt ready-to-dock molecules from .csv of SMILES strings.

        Parameters
        ----------
        which_ligands: string, ["gen", "train"] (Optional)
            Choose training or generated molecules to dock

        smiles_csv_path: Path or string (Optional)
            Path to CSV file containing SMILES strings enerated or train molecules)
            Note: overwrites `which_ligands` if specified
        """

        if which_ligands in ["gen", "train"]:
            self.which_ligands = which_ligands

        if smiles_csv_path == None:
            self.ligands_dir = self.data_path / self.paths["ligands_dir"]
            self.ligands_csv = (
                self.data_path / self.paths["ligands_dir"] / self.paths[f"ligands_{self.which_ligands}_csv"]
            )
            self.ligands_subdir = self.ligands_dir / self.paths[f"ligands_{self.which_ligands}_dir"
            os.makedirs(self.ligands_subdir, exists_ok=True)
        else:
            self.ligands_csv = Path(smiles_csv_path)

        self._assign_uuid()
        self._get_pdbqt()

    def _assign_uuid(self):
        """
        Appends a column of UUIDs to .csv of SMILES if the column isn't already there
        """

        ligands = pd.read_csv(ligands_csv)
        if "uuid" not in ligands:  # if UUIDs are not assigned yet
            ligands["uuid"] = ligands["SMILES"].apply(hash_smiles)

        ligands.to_csv(ligands_csv)

    def _get_pdbqt(self):
        """
        Generates PDBQT files from SMILES strings using open babel.
        """

        _ligands = pd.read_csv(ligands_csv)
        valid = [molecule_is_valid(molec) for molec in _ligands["SMILES"]]
        ligands = _ligands[valid].reset_index(drop=True)
        os.makedirs(self.ligands_subdir, exist_ok=True)

        for i in ligands.index:
            print(
                f"converting molecule {ligands['uuid'][i]} ({i + 1}/{len(ligands)})"
            )

            # write temporary SMILES
            molec = ligands["SMILES"][i]
            uuid = ligands["uuid"][i]
            temp_file = self.ligands_subdir / f"temp_{uuid}.smi"  # "valid_sample_1e5"
            with open(temp_file, "w") as temp:
                temp.write(molec)

            # TODO: openbabel conversion settings, use FullConvert might be faster
            conv = openbabel.OBConversion()
            conv.SetInAndOutFormats("smi", "pdbqt")
            conv.AddOption("gen3d", conv.GENOPTIONS)
            conv.AddOption("h", conv.GENOPTIONS)

            # convert temporary SMILES to PDBQT format
            pdbqt_file = self.ligands_subdir / f"{uuid}.pdbqt"
            conv.OpenInAndOutFiles(str(temp_file), str(pdbqt_file))
            conv.Convert()
            os.remove(temp_file)
        

    def dock(self, targets=None, subsample=1, dock_from_logs=None, rewrite=False):
        """
        Performs docking for the prepared molecules on specified targets.

        Parameters
        ----------
        targets: list (Optional)
            List of molecules to dock on. If not specified, then all will be docked

        subsample: float (Optional)
            Percentage of molecules to dock (0.0 to 1.0)

        dock_from_logs: string (Optional)
            Specifies path to a log file from which to take corresponding molecules for docking

        rewrite: bool (Optional)
            Re-dock molecules that have already been docked
        """

        if not targets:  # dock all
            targets_to_dock = list(self.targets["off_targets"].keys()) + list(
                self.targets["variants"].keys()
            )
        else:
            targets_to_dock = targets

        for trgt in targets_to_dock:
            trgt_name = trgt.lower()
            cur_out_dir = Path(
                self.data_path
                / self.paths["out_dir"]
                / self.paths[f"out_{self.which_ligands}_dir"]
                / trgt_name
            )
            out_pdbqt_dir = cur_out_dir / self.paths["out_pdbqt_dir"]
            out_logs_dir = cur_out_dir / self.paths["out_logs_dir"]
            os.makedirs(out_path, exist_ok=out_pdbqt_dir)
            os.makedirs(logs_path, exist_ok=out_logs_dir)

            ligands_to_dock = list(  # filter out already docked molecules
                set(os.listdir(self.ligands_gen_dir)) - set(os.listdir(out_pdbqt_dir))
            )

            if os.path.isfile(
                str(dock_from_logs)
            ):  # if choose molecules from logs of another docking experiment
                ligands_from_logs = pd.read_csv(molecs_from_logs)["uuid"].values
                ligand_files = []
                for lig in ligands_to_dock:
                    if change_file_ext(lig) in ligands_from_logs:
                        ligand_files.append(lig)

            elif subsample != False or subsample < 1:
                subsample_count = int(subsample * len(ligands_to_dock))
                print(f"subsampling {subsample_count} molecules to dock")
                ligand_files = sample(ligands_to_dock, subsample_count)

            else:
                ligand_files = ligands_to_dock

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

    def generate_logs_table(self, targets):
        pass
