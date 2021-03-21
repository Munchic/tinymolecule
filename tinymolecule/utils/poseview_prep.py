import os
from pathlib import Path

from openbabel import openbabel


DOCK_DIR = Path("/Users/Munchic/Developer/Capstone/tinymolecule/data/pdb")
PROT = "ccr5_valid_45e"
POSEVIEW_DIR = Path("/Users/Munchic/Developer/Capstone/tinymolecule/data/poseview")


def get_sdf_from_pdbqt(molec_path=DOCK_DIR / PROT, save_dir=POSEVIEW_DIR, mol_uuids=[]):
    os.makedirs(save_dir / PROT, exist_ok=True)
    save_folder = save_dir / PROT

    for i, mol_uuid in enumerate(mol_uuids):
        print(f"converting molecule {mol_uuid} ({i + 1}/{len(mol_uuids)})")

        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats("pdbqt", "sdf")

        pdbqt_file = molec_path / (mol_uuid + ".pdbqt")
        sdf_file = save_folder / (mol_uuid + "-hello.sdf")

        conv.OpenInAndOutFiles(str(pdbqt_file), str(sdf_file))
        conv.Convert()


get_sdf_from_pdbqt(
    mol_uuids=["ec0fdc31", "e8244ede", "82301e56", "802f05ac", "1f3586a9"]
)
