import os
from copy import deepcopy
from pathlib import Path
from typing import Optional, Union, List, Tuple

import yaml


class TinyTrain:
    """
    Part of the tinymolecule pipeline assisting in training generative models.
    This class implements data processing and filtering on downloaded activity data.
    """

    _DEFAULT_CONFIG = Path(__file__).parent / "config" / "default_config.yaml"

    def __init__(self, config=None):
        self.config = yaml.load(self._DEFAULT_CONFIG) if not config else config

    def get_train_data(self, smiles_csv: Union[string, Path]):
        pass
