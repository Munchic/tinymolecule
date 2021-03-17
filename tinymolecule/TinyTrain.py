import os
from copy import deepcopy
from pathlib import Path
from typing import Optional, Union, List, Tuple
import yaml

import torch
import pandas as pd
import moses


class TinyTrain:
    """
    Part of the tinymolecule pipeline assisting in training generative models.
    This class implements data processing and filtering on downloaded activity data.
    """

    _DEFAULT_CONFIG = Path(__file__).parent / "config" / "default_config.yaml"
    _ROOT_DIR = Path(__file__).parent.parent

    def __init__(self, config=None):
        with open(self._DEFAULT_CONFIG) as config_file:
            self.config = (
                yaml.load(config_file, Loader=yaml.Loader) if not config else config
            )

    def get_assays_data(self, assays_csv: Union[str, Path]):
        """
        Processes and filters assay data per config specification.

        Parameters
        ----------
        assays_csv: Path or string
            Path to a .csv containing SMILES strings
        """
        assays = pd.read_csv(assays_csv)
        self.paths = self.config["data_paths"]
        self.data_path = self._ROOT_DIR / self.paths["data_dir"]
        self.train_data_params = self.config["train_data_params"]

        self.ligands_csv = (
            self.data_path
            / self.paths["ligands_dir"]
            / self.paths[f"ligands_train_csv"]
        )

        self.assay_type = self.train_data_params["assay_type"]
        self.standard_units = self.train_data_params["standard_units"]
        self.smiles_type = self.train_data_params["smiles_type"]

        assays = self._filter_activity(assays)
        assays = self._filter_valid(assays)
        assays = self._rm_duplicates(assays)
        assays = self._rm_molecules(assays)
        assays.reset_index(inplace=True, drop=True)

        self.assays = assays
        print("assays data prepared, you can access at .assays")

    def write_train_csv(self):
        """
        Outputs a CSV containing just SMILES for training generative models
        """
        smiles_only = pd.DataFrame({"SMILES": list(self.assays[self.smiles_type])})
        smiles_only.to_csv(self.ligands_csv)

    def _filter_activity(self, assays):
        assays = assays[assays["type"] == self.assay_type]
        assays = assays[assays["standard_units"] == self.standard_units]

        return assays

    def _filter_valid(self, assays):
        assays = assays[assays[self.smiles_type].apply(lambda x: isinstance(x, str))]

        return assays

    def _rm_duplicates(self, assays):
        assays = assays.drop_duplicates()

        return assays

    def _rm_molecules(self, assays):
        smiles_to_rm = self.train_data_params["molecules_to_remove"].values()
        assays = assays[~assays[self.smiles_type].isin(smiles_to_rm)]

        return assays
