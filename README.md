# About
This method estimates the structural susceptibility of protein regions to proteolysis based on the known (from PDB database) or model (from AlphaFoldDB database) three-dimensional structure of a protein. The method calculates a numerical score, ranging from 0 to 1, for each peptide bond, representing its structural susceptibility to proteolysis. A score of 1 indicates a high structural susceptibility to proteolysis, while a score of 0 indicates a low structural susceptibility to proteolysis.

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

To calculate scores for a specific chain of the known experimental 3D structure (e.g. the A chain of 4GAW), execute the following command, which includes automatic downloading of the 3D structure from the PDB server (`-s`), selecting the chain (`-c`), running the locally installed DSSP (`-l`), and visualizing the results (`-v`). The known structure ID must be four-symbol, according to entry ID from PDB database. Case-insensitive form is permitted. The name of the protein chain is case-sensitive. Please, check the protein chain you are interested in before calculating proteolytic scores.

```python run.py -s 4GAW -c A -l -v```
```python run.py -s 1abe -c A -l -v```

To apply the method for locally saved PDB-file of an experimental structure, also use `-s` option.

```python run.py -s /dir1/dir2/my_experimental_protein.pdb -c A -l -v```

To calculate scores for a model 3D structure from AlphaFoldDB (e.g. P00257), execute the following command, which includes automatic downloading of the 3D structure from the AlphaFoldDB server (`-m`), running the locally installed DSSP (`-l`), and visualizing the results (`-v`). The model structure ID is a UniProt ID. Selecting of the protein chain is not required in this situation.

```python run.py -m P00257 -l -v```

To apply the method for locally saved PDB-file of a model structure, also use `-m` option.

```python run.py -m /dir1/dir2/my_model_protein.pdb -l -v```

For batch execution of the same command on multiple structures (including mix set - experimental and model structures, ID of structures and locally saved PDB-files), you can submit a text file with a list of structures and chains, with each entry in two columns (the ID and the protein chain) per line (It's recommended to separate values by comma). Specify '-' for AlphaFold structure chain value. Importantly, in this option the name of your model structure ID must match a template as "AF-{UniProt_ID}-F1", and the name of your model structure PDB-file must start with "AF-{your_name}.pdb" for correct processing:

```python run.py -b "/dir1/dir2/my_structures.txt" -l -v```

Please, see `input_examples` directory for better understanding.

The calculation results are accessible in the `results` directory within a folder named after the corresponding structure ID. The calculated scores in a concise format can be obtained in the file `<Structure_ID>.predictions.txt`.

# User manual

All calculations are performed using the main script `run.py`. To obtain a full list of options, you can execute the `python run.py --help` command.

* `-s` - the four-symbol structure ID from PDB database (case-insensitive), or a full path to the known structure file with '.pdb' extension. It must be used in conjunction with the '-c' parameter. For example, 1) -s '4gaw' -c 'A' 2) -s '/dir1/my_experimental_structure.pdb' -c 'A'. Also see directory 'input_examples' for a better understanding.

* `-m` - the UniProt ID of a protein for downloading model structure from AlphaFoldDB database, or a full path to the model structure file with '.pdb' extension. For example, 1) -m 'P00257' 2) -m '/dir1/dir2/my_model_structure.pdb'. Also see directory 'input_examples' for a better understanding.

* `-b` - a batch option to generate predictions for a set of structures. A full path to the TXT-file with a table where the first column is the structure ID (The name of your AlphaFold structure file must be necessarily started with 'AF-' and the name of your AlphaFold structure ID must match a template as "AF-{UniProt_ID}_F1"), and the second is the name of the chain (You need specify '-' for AlphaFold structure). The value separator should be a comma. For example, 1) -b '/dir1/dir2/structures.txt'. Also see directory 'input_examples' for a better understanding.

* `-c` - the name of the chain (case-sensitive), only used with '-s' parameter for PDB ID or PDB structure file. For example, 'A'."

* `-l` (optional) - run DSSP program locally. It is strongly recommended to install DSSP locally due to potential issues with the connection to the DSSP server.

* `-v` (optional) - generate a visual representation of the 3D structure in Chimera with color-mapped scores. The colors reflecting the scores of structural susceptibility to proteolysis are distributed from blue (low susceptibility) through red (moderate susceptibility) to yellow (high susceptibility). The visual representation file in Chimera format is placed in the `scores` folder.

After completing the calculations, the main script creates a directory named `results` if it does not already exist, and places the final and intermediate output files under the folder with the structure ID name:

* `<Structure_ID>.predictions.txt` - the calculated scores of structural susceptibility to proteolysis assigned to the P1 site (Schechter-Berger notation) in a file with a concise format.

* `scores` - folder containing detailed output of the method, including the predicted structural susceptibility scores and all associated structural features. If visualization option was enabled, this folder will also contain a Chimera session file.

* `features` - folder containing structural features before and after normalization.

* `alignments` - folder containing the results of auxiliary alignments for correct extraction of structural features.

* `preprocessing` - folder containing the DSSP-file and the results of processing of PDB-files (removing of water and ligand molecules).

Along with the `results` folder, the root directory also contains system folders `models` and `scripts`, which should not be removed or edited.
 

# Reporting Bugs and Feature Requests
Please use the [GitHub issue tracker](https://github.com/KazanovLab/ProteolysisStructuralPrediction/issues) to report bugs or suggest features.

# Citing
Matveev, E.V.; Safronov, V.V.; Ponomarev, G.V.; Kazanov, M.D. Predicting Structural Susceptibility of Proteins to Proteolytic Processing. Int. J. Mol. Sci. 2023, 24, 10761. https://doi.org/10.3390/ijms241310761
