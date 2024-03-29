# About
This method estimates the structural susceptibility of protein regions to proteolysis based on the known three-dimensional structure of a protein. The method calculates a numerical score, ranging from 0 to 1, for each peptide bond, representing its structural susceptibility to proteolysis. A score of 1 indicates a high structural susceptibility to proteolysis, while a score of 0 indicates a low structural susceptibility to proteolysis.

# Dependencies
For the correct work of the method, you need to install the following software:
* [Clustal Omega](http://www.clustal.org/omega/)
* [UCSF Chimera](https://github.com/insilichem/pychimera/blob/master/docs/install.rst)
* [DSSP](https://github.com/PDB-REDO/dssp) (optional, but strongly recommended)

and the following python libraries:

* Biopython
* pandas
* numpy
* requests

# Quick start

To calculate scores for a specific chain of the 3D structure (e.g. the A chain of 4GAW), execute the following command, which includes automatic downloading of the 3D structure from the PDB server (`-i`), selecting the chain (`-c`), running the locally installed DSSP (`-l`), and visualizing the results (`-v`):

```python run.py -i 4GAW -c A -l -v```

For batch execution of the same command on multiple structures, you can submit a text file with a list of structures and chains, with each entry in two columns (PDB ID and Chain ID) per line:

```python run.py -i "PDBID_ChainID_list.txt" -l -v```

To apply the method for locally saved PDB file, you should use `-f` option:

```python run.py -f 4GAW.pdb -c A -l -v```

The batch option is also available for locally saved PDB files listed in a two-column text file, with the first column representing the path to the PDB file and the second column representing the Chain ID.

```python run.py -f "PDBpath_ChainID_list.txt" -l -v```

The calculation results are accessible in the `results` directory within a folder named after the corresponding PDB ID. The calculated scores in a concise format can be obtained in the file `<PDBID_ChainID>_predictions.txt`.

# User manual

All calculations are performed using the main script `run.py`. To obtain a full list of options, you can execute the `python run.py --help` command.

* `-i` - PDB ID of the input 3D structure for downloading from PDB server or the path to the batch file. The batch file is a text file with a two-column format (.txt extension), where the first column contains the PDB ID and the second column contains the Chain ID.

* `-f` - path to the local PDB file or the path to the batch file. The batch file in this case also consists of two columns: the first column contains the path to the PDB file, and the second column contains the Chain ID.

* `-c` - chain ID

* `-l` (optional) - run DSSP locally. It is strongly recommended to install DSSP locally due to potential issues with the connection to the DSSP server.

* `-v` (optional) - generate a visual representation of the 3D structure in Chimera with color-mapped scores. The colors reflecting the scores of structural susceptibility to proteolysis are distributed from blue (low susceptibility) through red (moderate susceptibility) to yellow (high susceptibility). The visual representation file in Chimera format is placed in the `scores` folder.

After completing the calculations, the main script creates a directory named `results` if it does not already exist, and places the final and intermediate output files under the folder with the PDB ID name:

* `<PDBID_ChainID>_predictions.txt` - the calculated scores of structural susceptibility to proteolysis assigned to the P1 site (Schechter-Berger notation) in a file with a concise format.

* `scores` - folder containing detailed output of the method, including the predicted structural susceptibility scores and all associated structural features. If visualization option was enabled, this folder will also contain a Chimera session file.

* `features` - folder containing structural features before and after normalization.

Along with the `results` folder, the root directory also contains system folders `models` and `scripts`, which should not be removed or edited.
 

# Reporting Bugs and Feature Requests
Please use the [GitHub issue tracker](https://github.com/KazanovLab/ProteolysisStructuralPrediction/issues) to report bugs or suggest features.

# Citing
Matveev, E.V.; Safronov, V.V.; Ponomarev, G.V.; Kazanov, M.D. Predicting Structural Susceptibility of Proteins to Proteolytic Processing. Int. J. Mol. Sci. 2023, 24, 10761. https://doi.org/10.3390/ijms241310761
