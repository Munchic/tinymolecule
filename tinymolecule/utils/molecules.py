import numpy as np

import moses
from moses.metrics.utils import (
    get_mol,
    mol_passes_filters,
    logP,
    QED,
    SA,
    weight,
    get_n_rings,
)


def molecular_properties(molec: str):
    mol_obj = get_mol(molec)
    return {
        "MCF_PAINS_pass": mol_passes_filters(mol_obj),
        "logP": logP(mol_obj),
        "QED": QED(mol_obj),
        "SA": SA(mol_obj),
        "weight": weight(mol_obj),
        "n_rings": get_n_rings(mol_obj),
    }


def molecule_is_valid(molec: str):
    return get_mol(molec) not in [None, np.nan]