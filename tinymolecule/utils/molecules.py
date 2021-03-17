import numpy as np

import moses


def molecule_is_valid(molec: str):
    return moses.utils.get_mol(molec) not in [None, np.nan]