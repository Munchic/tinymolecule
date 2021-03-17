import moses


def molecule_is_valid(molec: string):
    return moses.utils.get_mol(molec) not in [None, np.nan]