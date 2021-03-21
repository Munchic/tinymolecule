# *tinymolecule*: in silico experiments on off-target binding and adaptability to genetic variants of small molecules generated with deep learning
## Installation
1. Install [AutoDock Vina CLI](http://vina.scripps.edu/manual.html#installation) for performing docking experiments
2. Create a `conda` environment and inside it:
    ```bash
    $ conda create --name tinymolecule python=3.6
    $ conda activate tinymolecule
3. Install [RDKit](https://www.rdkit.org/docs/Install.html) and [MOSES](https://github.com/molecularsets/moses) inside the `conda` environment for using MOSES metrics and performing molecule validity checks
4. Install [openbabel](https://pypi.org/project/openbabel/) in the `conda` environment for molecule format conversions
    ```bash
    $ pip install openbabel
    ```

## Quick start
There are three main classes that you interact with in `tinymolecule`: `TinyTrain()`, `TinyDock()`, and `TinyAnalyze()`. If you just want to dock your molecules, skip to section II. 

### I. Preparing training dataset for generative models


### II. Processing ligands and high-throughput docking


### III. Ranking docked ligands and downstream analyses