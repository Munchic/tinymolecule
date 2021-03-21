# *tinymolecule*: in silico experiments on off-target binding and adaptability to genetic variants of small molecules generated with deep learning
## Installation
1. Clone *tinymolecule* repository:
    ```bash
    $ git clone https://github.com/Munchic/tinymolecule.git
    ```
2. Install [AutoDock Vina CLI](http://vina.scripps.edu/manual.html#installation) for performing docking experiments
3. Create a `conda` environment and inside it:
    ```bash
    $ conda create --name tinymolecule python=3.6
    $ conda activate tinymolecule
4. Install [RDKit](https://www.rdkit.org/docs/Install.html) and [MOSES](https://github.com/molecularsets/moses) inside the `conda` environment for using MOSES metrics and performing molecule validity checks
5. Install [openbabel](https://pypi.org/project/openbabel/) in the `conda` environment for molecule format conversions
    ```bash
    $ pip install openbabel
    ```
## Quick start
### I. Playground
1. Download the whole dataset containing processed targets, ligands, as well as docking results from [here](https://drive.google.com/file/d/1uht6BvRad-liuu7nqa4p_Jnya4opA8QE/view?usp=sharing) and unpack directly into the repository. Hint: `data` folder should be on the same level as `analysis`
2. `ligands` folder contains `.csv` summaries of what training and generated molecules we have, and the folders store the dock-able `.pdbqt` files (more info in Section II.2)
3. `targets` folder contains prepared 3D structures of protein targets
4. `out` folder contains docked molecules
5. Follow Section I.1-I.2
6. There are a few operations we can do with existing data, for example, dock a small set of molecules from the `ligands` folder:
    ```python
    td = tinymolecule.TinyDock()
    td.dock(targets=["WT_CCR5"], subsample=0.005)
    ```
    We will dock 0.5% of the ligands onto wild-type CCR5
7. You will see that the folders `data/out/gen/wt_ccr5/pdbqt/` and `data/out/gen/wt_ccr5/logs` start updating with new docking results
8. To view summaries of logs, you can open `data/out/gen/wt_ccr5/summary.csv`
9. Let's try to prioritize some molecules and see what they look like:
    ```python
    ta = tinymolecule.TinyAnalyze()
    ta.prioritize()  # sort molecules based on objective function
    ta.priority_df.head()  # view top molecules 
    ```

    ```python
    smiles = ta.get_smiles_from_uuid("<PUT_UUID_HERE>")  # fetch corresponding SMILES

    from moses.utils import get_mol
    get_mol(smiles)
    ```
    Ta-da, here you have a generated molecule!
10. For more detailed guidance, follow the Instructions


### II. Analysis notebooks
There are a few analysis notebooks used to produce figures in the paper located in `analysis/`
1. `baff_gen_vs_train.ipynb` compares generated molecules to training
2. `baff_off_target.ipynb` compares binding affinity on off-targets vs CCR5
3. `baff_gen_vs_train.ipynb` compares binding affinity on CCR5 variants vs wild-type CCR5
4. `baff_best.ipynb` prioritizes the top five molecules

## Instructions
There are three main classes that you interact with in *tinymolecule*: `TinyTrain()`, `TinyDock()`, and `TinyAnalyze()`. The central place to configure paths and parameters is `tinymolecule/config/default_config.yaml`. If you just want to dock your molecules, skip to Section II. 

### I. Preparing training dataset for generative models
1. Create a notebook inside the `analysis/` folder
2. Add the following cell to import *tinymolecule*:
    ```python
    import os
    import sys
    sys.path.append("..")  # add top folder to path

    import tinymolecule
    ```
3. Download assays metadata from ChEMBL (you might need to `pip install chembl_webresource_client` to install the ChEMBL web slient):
    ```python
    tinymolecule.utils.download_chembl_assays_metadata(TARGET_CHEMBL_ID, DIRECTORY_TO_SAVE, CSV_NAME)
    ```
    For example:
     ```python
    tinymolecule.utils.download_chembl_assays_metadata("CHEMBL274", "../data", "meta.csv")  # download assay data on CCR5
    ```
4. Open `default_config.yaml` and observe `# TRAINING DATA PROCESSING`
    1. Here you can configure molecules you want to remove
    2. Assay type and units
    3. Type of SMILES you want to retrieve
    4. (Optional) Modify the parameters to suit your use case
5. Back to the notebook, now we can run molecules processing with specified parameters to obtain a training set containing SMILES
    ```python
    tt = tinymolecule.TinyTrain()  # initialize
    tt.get_assays_data(PATH_TO_META_CSV)  # input the ChEMBL assays csv here
    tt.write_train_csv()  # create a training dataset
    ```
    Note that the training dataset will be stored at key `ligands_train_csv` parameter configured in `default_config.yaml`. This means you don't need to pass any parameters there, just make sure that the config is set correctly.
6. Your training dataset is ready at directory specified by the `ligands_train_csv` key
7. After you train your generative model and 

### II. Processing ligands and high-throughput docking
1. Do steps 1-2 in Section I if you haven't
2. Prepare your protein targets:
    1. Remove water and any ligands, add charged hydrogens
    2. Save as`.pdbqt` and name it `<protein_name>.pdbqt` similar to the format in `default_config.yaml`. Store this file in `data/targets/pdbqt/`
    3. Create `vina_config_<protein_name>.text` file containing search box coordinates and store it in `data/targets/vina_config/`
3. Prepare `.pdbqt` files from SMILES strings to dock
    ```python
    td = tinymolecule.TinyDock()
    td = tinymolecule.prepare_molecules(which_ligands="gen")
    ```
    `which_ligans` can either be `"train"` or `"gen"` depending on whether you're working with training data or generated molecules. For example, if you choose `"gen"`, when you process the molecules, they will be stored in where the key `ligands_gen_dir` specifies.

4. Let's dock the molecules! 
    ```python
    td.dock(targets=["TARGET_NAME"], subsample=1)
    ```
    Targets specifies which target from the config you would like to dock on. Make sure that the names match exactly to the key you specified in `default_config.py`. Parameter `subsample` chooses a random subsample to dock (that is, if you want to test if things are working).

    Example output:
    ```python
    >>> >>> DOCKING 6036 MOLECULES ON CCR2 <<< 
    >>> subsampling 6036 molecules to dock
    >>> docking molecule 2de25107 (1/6036) ❌ subprocess error
    >>> docking molecule edb0ba0a (2/6036) ✅ success
    >>> docking molecule ace426bc (3/6036) ✅ success
    >>> docking molecule af4a11ad (4/6036) ✅ success
    >>> docking molecule e2fedb66 (5/6036) ✅ success
    >>> docking molecule 46c9ee4a (6/6036) ✅ success
    >>> docking molecule 05208265 (7/6036)
    ```
5. The docked `.pdbqt` files will appearing under the directory specified by the key `out_pdbqt_dir` and the logs that include binding affinity should appear in `out_logs_dir`

### III. Ranking docked ligands and downstream analyses
1. Summarize the log files into a `summary.csv` for each target:
    ```
    ta = tinymolecule.TinyAnalyze()
    ta.summarize_logs()
    ```
    Summary files will appear at the path specified by `summary_csv`
2. Now, the juice of all, we can prioritize the molecules based on low binding affinity towards off-targets and high binding affinity towards on-targets:
    ```python
    ta.prioritize()
    ta.priority_df  # view dataframe of prioritized molecules
    ```
3. We can also look at the molecular properties from MOSES:
    ```python
    ta.get_molecular_properties(["ec0fdc31",
    "e8244ede"])
    ```
    The input is a list of UUIDs of the molecules.

    Sample output:
    ```
    {'ec0fdc31': {'MCF_PAINS_pass': True,
        'logP': 5.994200000000007,
        'QED': 0.24949291816662938,
        'SA': 2.996830627015026,
        'weight': 577.7910000000002,
        'n_rings': 4},
     'e8244ede': {'MCF_PAINS_pass': True,
        'logP': 5.825800000000007,
        'QED': 0.3397392620058908,
        'SA': 2.6949474504507673,
        'weight': 490.71300000000025,
        'n_rings': 4},
    }
    ```