import os
from copy import deepcopy
from pathlib import Path
from typing import Optional, Union, List, Tuple

import yaml


class TinyDock:
    """
    Core class of the tinymolecule pipeline to dock generated molecules.
    Stores docked molecules objects and their metrics. Implements analyses
    showcased in the paper.
    """

    _DEFAULT_CONFIG = Path(__file__).parent / "config" / "default_config.yaml"

    def __init__(self, config=None):
        self.config = yaml.load(self._DEFAULT_CONFIG) if not config else config

    def prepare_molecules(self, smiles_csv: Union[string, Path]):
        pass

    def dock(self, target=None, rewrite=False):
        '''
        
        '''

        if not target:  # dock all
            pass
        else:  # off-target or variants
            pass  # take from config    

    def 