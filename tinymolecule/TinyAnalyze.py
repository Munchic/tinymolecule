import os
from copy import deepcopy
from pathlib import Path
from typing import Optional, Union, List, Tuple
from functools import reduce

import yaml
import pandas as pd

from moses.metrics.utils import get_mol, logP, weight, get_n_rings, mol_passes_filters

from tinymolecule.utils.molecules import molecular_properties
from tinymolecule.utils.helper_fn import get_mean_baff_df


class TinyAnalyze:
    """
    This class implements analyses, small molecule prioritization based on
    objective function, and plotting.
    """

    _DEFAULT_CONFIG = Path(__file__).parent / "config" / "default_config.yaml"
    _ROOT_DIR = Path(__file__).parent.parent

    def __init__(self, config=None):
        with open(self._DEFAULT_CONFIG) as config_file:
            self.config = (
                yaml.load(config_file, Loader=yaml.Loader) if not config else config
            )
        self.paths = self.config["data_paths"]
        self.data_path = self._ROOT_DIR / self.paths["data_dir"]
        self.targets = self.config["targets"]

        self.assays_meta_csv = self.data_path / self.paths["assays_meta_csv"]

        self.ligands_train_dir = (
            self.data_path / self.paths["ligands_dir"] / self.paths["ligands_train_dir"]
        )
        self.ligands_gen_dir = (
            self.data_path / self.paths["ligands_dir"] / self.paths["ligands_gen_dir"]
        )
        self.ligands_train_csv = (
            self.data_path
            / self.paths["ligands_dir"]
            / self.paths[f"ligands_train_csv"]
        )
        self.ligands_gen_csv = (
            self.data_path / self.paths["ligands_dir"] / self.paths[f"ligands_gen_csv"]
        )

        self.out_train_subdir = (
            self.data_path / self.paths["out_dir"] / self.paths[f"out_train_dir"]
        )
        self.out_gen_subdir = (
            self.data_path / self.paths["out_dir"] / self.paths[f"out_gen_dir"]
        )

    def summarize_logs(self):  # TODO: fix paths
        """
        Prepares a summary file in each output directory, concatenating all the logs for that target.
        """

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

    def plot_binding_affinity_distribution(self, targets: list = None):
        """
        Generates an overlay of plots showing distributions of binding affinities across specified targets

        Parameters
        ----------
        targets: list of strings (Optional)
            Targets to show on the plot
        """

        if not targets:  # show all
            targets_to_show = list(self.targets["off_targets"].keys()) + list(
                self.targets["variants"].keys()
            )
        else:
            targets_to_show = targets

    def plot_binding_affinity_scatter(self, target_x: str, target_y: str):
        """
        Generates a scatter plot showing comparing binding affinities ligand-wise on two targets.
        Additionally, calculates Pearson's r and fits a linear regression for the relationship.

        Parameters
        ----------
        target_x: string
            Targets for which to display binding affinity on x-axis

        target_y: string
            Targets for which to display binding affinity on y-axis
        """

        pass

    def get_mol_props_and_baff(self, which_ligands="gen", target="WT_CCR5"):
        """
        Calculates a dataframe of molecular properties in correspondence with binding affinities.

        | uuid | logP | weight | n_rings | MCF | affin_kcal_mol-1_avg | affin_kcal_mol-1_1 | ... | affin_kcal_mol-1_n |
        """

        target_name = target.lower()
        summary_path = self.out_gen_subdir / target_name / self.paths["summary_csv"]
        baff_df = pd.read_csv(summary_path)
        ligands_gen_df = pd.read_csv(self.ligands_gen_csv)

        # merge summary table and gen table
        mol_props_and_baff_df = pd.merge(baff_df, ligands_gen_df, how="left", on="uuid")

        # use smiles to calculate properties
        mol_objs = [get_mol(smi) for smi in mol_props_and_baff_df["SMILES"]]
        mol_props_and_baff_df["logP"] = [logP(mol) for mol in mol_objs]
        mol_props_and_baff_df["weight"] = [weight(mol) for mol in mol_objs]
        mol_props_and_baff_df["n_rings"] = [get_n_rings(mol) for mol in mol_objs]
        mol_props_and_baff_df["MCF_PAINS"] = [
            mol_passes_filters(mol) for mol in mol_objs
        ]

        return mol_props_and_baff_df

    def prioritize(self, which_ligands="gen"):
        """
        Returns a dataframe of small molecules prioritized using the objective function

        | uuid | variants_mean | off_target_max | objective |
        """

        logs_variants = {}
        for target in self.targets["variants"]:
            target_name = target.lower()
            cur_summary = pd.read_csv(self.out_gen_subdir / target_name / "summary.csv")
            logs_variants[target_name] = get_mean_baff_df(
                cur_summary, prefix=f"{target_name}_"
            )

        logs_off_targets = {}
        for target in self.targets["off_targets"]:
            target_name = target.lower()
            cur_summary = pd.read_csv(self.out_gen_subdir / target_name / "summary.csv")
            logs_off_targets[target_name] = get_mean_baff_df(
                cur_summary, prefix=f"{target_name}_"
            )

        variants_merged = reduce(
            lambda left, right: pd.merge(left, right, on="uuid"), logs_variants.values()
        )
        off_targets_merged = reduce(
            lambda left, right: pd.merge(left, right, on="uuid"),
            logs_off_targets.values(),
        )

        # the norms
        variants_merged["variants_mean"] = variants_merged.mean(axis=1)
        off_targets_merged["off_targets_max"] = off_targets_merged.max(axis=1)

        cross_baff = pd.merge(off_targets_merged, variants_merged, on="uuid")
        cross_baff["objective"] = (
            cross_baff["variants_mean"] - cross_baff["off_targets_max"]
        )
        cross_baff.sort_values(by="objective", ascending=False, inplace=True)

        self.priority_df = cross_baff[
            ["uuid", "variants_mean", "off_targets_max", "objective"]
        ]
        print("molecules prioritized, you can access at TinyAnalyze.priority_df")

    def get_molecular_properties(self, molecules_uuid: list):
        """
        For a given list of molecular UUIDs, returns a dictionary of molecular properties.
        """

        props_dict = {}
        ligands_gen_df = pd.read_csv(self.ligands_gen_csv)

        # TODO: add uuid exists in .csv check
        for m_uuid in molecules_uuid:
            smiles = ligands_gen_df[ligands_gen_df["uuid"] == m_uuid][
                "SMILES"
            ].reset_index(drop=True)
            props = molecular_properties(smiles[0])
            props_dict[m_uuid] = props

        return props_dict

    def get_smiles_from_uuid(self, uuid: str):
        ligands_gen_df = pd.read_csv(self.ligands_gen_csv)
        smiles = ligands_gen_df[ligands_gen_df["uuid"] == uuid]["SMILES"].reset_index(
            drop=True
        )

        return smiles[0]
