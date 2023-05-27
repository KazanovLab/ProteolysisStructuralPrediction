# About
This method estimates the structural susceptibility of protein regions to proteolysis based on the known three-dimensional structure of a protein. The method calculates a numerical score, ranging from 0 to 1, for each peptide bond, representing its structural susceptibility to proteolysis. A score of 1 indicates a high structural susceptibility to proteolysis, while a score of 0 indicates a low structural susceptibility to proteolysis.

# Dependencies
For the correct work of the method, you need to install the following software:
* [Clustal Omega](http://www.clustal.org/omega/)
* [DSSP](https://github.com/PDB-REDO/dssp) 
* [UCSF Chimera](https://github.com/insilichem/pychimera/blob/master/docs/install.rst)

and the following python libraries:

* Biopython
* pandas
* numpy
* requests

# Quick start

To calculate scores for a specific chain of the 3D structure (e.g. the A chain of 4GAW), execute the following command, which includes automatic downloading of the 3D structure from the PDB server (`-i`), running the locally installed DSSP (`-l`), and visualizing the results (`-v`):

```python run.py -i "4GAW_A" -l -v```

For batch execution of the same command on multiple structures, you can submit a text file with a list of structures and chains, with one entry per line:

```python run.py -i "ID_list.txt" -l -v```

# User manual

# Reporting Bugs and Feature Requests
Please use the [GitHub issue tracker](https://github.com/KazanovLab/ProteolysisStructuralPrediction/issues) to report bugs or suggest features.

# Citing
_Evgenii V. Matveev, Vyacheslav V. Safronov, Gennady V. Ponomarev and Marat D. Kazanov. "Predicting structural susceptibility of proteins to proteolytic processing". (submitted)_