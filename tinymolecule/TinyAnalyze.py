import os
from copy import deepcopy
from pathlib import Path
from typing import Optional, Union, List, Tuple

import yaml


class TinyAnalyze:
    """
    Core class of the tinymolecule pipeline to dock generated molecules.
    Stores docked molecules objects and their metrics. Implements analyses
    showcased in the paper.
    """

    _DEFAULT_CONFIG = Path(__file__).parent / "config" / "default_config.yaml"

    def __init__(self, config=None):
        self.config = yaml.load(self._DEFAULT_CONFIG) if not config else config
        self.paths = config["data_paths"]
        self.data_path = Path(config["data_paths"]["data"])
        self.targets = config["targets"]

    def prepare_molecules(
        self, which_data: "gen", smiles_csv_path: Union[str, Path] = None
    ):
        """
        Creates pdbqt ready-to-dock molecules from .csv of SMILES strings.

        Parameters
        ----------
        which_data: string, ["gen", "train", "all"] (Optional)
            Choose training or generated molecules to dock

        smiles_csv_path: Path or string (Optional)
            Path to CSV file containing SMILES strings enerated or train molecules)
            Note: overwrites `which_data` if specified
        """

        if smiles_csv_path == None:
            self.gen_csv = self.data_path / self.paths[f"ligands_{dataset}_csv"]
        else:
            self.gen_csv = Path(smiles_csv_path)

    def dock(self, target=None, rewrite=False):
        """"""

        if not target:  # dock all
            pass
        else:  # off-target or variants
            pass  # take from config
