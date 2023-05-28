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

To calculate scores for a specific chain of the 3D structure (e.g. the A chain of 4GAW), execute the following command, which includes automatic downloading of the 3D structure from the PDB server (`-i`), selecting the chain (`-c`), running the locally installed DSSP (`-l`), and visualizing the results (`-v`):

```python run.py -i 4GAW -c A -l -v```

For batch execution of the same command on multiple structures, you can submit a text file with a list of structures and chains, with each entry in two columns (PDB ID and Chain ID) per line:

```python run.py -i "PDBID_ChainID_list.txt" -l -v```

To apply the method for locally saved PDF file or multiple PDB files in folder, you should use `-f` option:

```python run.py -f 4GAW.pdb -c A -l -v```

The batch option is also available for locally saved PDB files listed in a two-column text file, with the first column representing the path to the PDB file and the second column representing the Chain ID.

```python run.py -i "PDBpath_ChainID_list.txt" -l -v```

The calculation results are accessible in the `results` directory within a folder named after the corresponding PDB ID. The calculated scores in a concise format can be obtained in the file `<PDB_ID>.scores.txt`.

# User manual

All calculations are performed using the main script `run.py`. To obtain a full list of options, you can execute the `python run.py --help` command.

* `-i` - PDB ID of the input 3D structure for downloading from PDB server or the path to the batch file. The batch file is the two-columns text file (`.txt` extention) with 
* `-f` - path to the local PDB file or the path to the batch file.

# Reporting Bugs and Feature Requests
Please use the [GitHub issue tracker](https://github.com/KazanovLab/ProteolysisStructuralPrediction/issues) to report bugs or suggest features.

# Citing
_Evgenii V. Matveev, Vyacheslav V. Safronov, Gennady V. Ponomarev and Marat D. Kazanov. "Predicting structural susceptibility of proteins to proteolytic processing". (submitted)_