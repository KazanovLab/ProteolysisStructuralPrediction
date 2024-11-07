""" Import libraries """
from Bio.PDB import *
import requests
import pandas as pd

import os
import sys
import time
import shutil
import argparse

''' Terminal '''
parser = argparse.ArgumentParser()
parser.add_argument("-s", "--structure", help="The four-letter structure ID from PDB database (case-insensitive), or a full path to the structure file with '.pdb' extension. It must be used in conjunction with the '-c' parameter. For example, 1) -s '4GAW' -c 'A' 2) -s '/dir1/my_experiment_structure.pdb' -c 'A'. Also see directory 'input_examples' for a better understanding.", type=str)
parser.add_argument("-m", "--model", help="The UniProt ID of a protein for downloading model structure from AlphaFoldDB database, or a full path to the model structure file with '.pdb' extension. For example, 1) -m 'P00257' 2) -m '/dir1/dir2/my_model_structure.pdb'. Also see directory 'input_examples' for a better understanding.", type=str)
parser.add_argument("-b", "--batch", help="A batch option to generate predictions for a set of structures. A full path to the TXT-file with a table where the first column is the structure ID (The name of your AlphaFold structure file must be necessarily started with 'AF-'), and the second is the name of the chain (You need specify '-' for AlphaFold structure). The value separator should be a comma. For example, 1) -b '/dir1/dir2/structures.txt'. Also see directory 'input_examples' for a better understanding.", type=str)
parser.add_argument("-c", "--chain", help="A name of the chain, only for PDB ID or PDB file. For example, 'A'.", type=str)
args = parser.parse_args()

''' Set paths '''
main_path = os.getcwd()
StructureSet_path = os.path.join(main_path, "results")
if not os.path.exists(StructureSet_path):
    os.mkdir(StructureSet_path)

''' Logging '''    
log_file = os.path.join(StructureSet_path, "report.log")

''' List of structures '''
# PDB #  
if (not args.structure is None) and (args.model is None) and (args.batch is None):
    structure_list = [[args.structure, args.chain]]
# AlphaFold #
elif (args.structure is None) and (not args.model is None) and (args.batch is None):
    structure_list = [[args.model, '-']]
# Batch #
elif (args.structure is None) and (args.model is None) and (not args.batch is None):
    with open(args.batch, 'r') as file:
        line = file.readline()
        if ',' in line:
            sep_value = ','
        elif ';' in line:
            sep_value = ';'
        elif '\t' in line:
            sep_value = '\t'
        elif (' ' in line) and (len(line.split(' ')) == 2):
            sep_value = ' '
        else:
            print("ERROR: Unrecognized separator value in table!\nPlease, use ',', ';', '\\t', ' ' as separator value.\n")
            sys.exit()
    input_data = pd.read_csv(args.batch, names=["structure", "chain"], sep=sep_value, dtype={"structure":str, "chain":str})
    structure_list = input_data.values
else:
    print("ERROR: Possibly, your input parameters are uncorrected! Check with templates:\n1: -s 4GAW -c A\n2: -s /dir1/my_experimental_structure.pdb -c A\n3: -m P00257\n4: -m /dir1/dir2/my_model_structure.pdb\n5: -b /dir1/dir2/my_structures.txt\n")
    sys.exit()
     
def download_PDB(structure_id):
	response = requests.get(f"https://files.rcsb.org/download/{structure_id}.pdb")
	return response.text

def download_AF(structure_id):
	response = requests.get(f"https://alphafold.ebi.ac.uk/files/{structure_id}-model_v4.pdb")
	return response.text

def extract_chain(path_to_read, structure_with_chain, path_to_save):
    class ChainSelect(Select):

        def __init__(self, needed_chain):
            self.needed_chain = needed_chain

        def accept_chain(self, chain):
            if chain.id == self.needed_chain:
                return 1
            else:
                return 0

    structure_id = '_'.join(structure_with_chain.split('_')[:-1])
    structure_chain = structure_with_chain.split('_')[-1]
    
    parser = PDBParser()
    io = PDBIO()
    structure_io = parser.get_structure(structure_id, f"{path_to_read}/{structure_id}.pdb")
    if len(structure_io) > 1:
        print("CHECK (for developers)!")
        sys.exit()
    else:
        io.set_structure(structure_io)
        io.save(f"{path_to_save}/{structure_with_chain}.pdb", ChainSelect(structure_chain))

def main():
    uncorrected_structure_files = []
    
    for num, (structure, chain) in enumerate(structure_list):
        #print(f"{num + 1} --- {structure}")
        
        if os.path.isfile(structure):
            if args.model is None:
                basename_structure_file = os.path.basename(structure)
                structure_id = basename_structure_file.split('.pdb')[0]
            else:
                basename_structure_file = os.path.basename(structure)
                if basename_structure_file.split('.pdb')[0].startswith('AF-') and basename_structure_file.split('.pdb')[0].endswith('-F1'):
                    structure_id = basename_structure_file.split('.pdb')[0]
                else:
                    structure_id = f"AF-{basename_structure_file.split('.pdb')[0]}-F1"
        else:
            if args.model is None:
                structure_id = structure
            else:
                structure_id = f"AF-{structure}-F1"
        
        structure_path = os.path.join(StructureSet_path, structure_id)
        if not os.path.exists(structure_path):
            os.mkdir(structure_path)
        
        preprocessing_path = os.path.join(structure_path, 'preprocessing')
        if not os.path.exists(preprocessing_path):
            os.mkdir(preprocessing_path)
             
        ''' Download structure file if required'''
        if os.path.isfile(structure):
            shutil.copyfile(structure, os.path.join(structure_path, basename_structure_file))  
        else:
            structure_file = os.path.join(structure_path, f"{structure_id}.pdb")
            if not os.path.exists(structure_file):
                structure_data = download_PDB(structure_id) if "AF-" not in structure_id else download_AF(structure_id)
                if "ATOM" not in structure_data:
                    uncorrected_structure_files.append(f"{structure}")
                    shutil.rmtree(structure_path)
                    continue
                else:
                    list_to_write = []
                    for stroka in structure_data.strip('\n').split('\n'):
                        list_to_write.append(stroka)
                        if stroka.split(' ')[0] == "ENDMDL":
                            break
                    with open(structure_file, "w") as file:
                        file.write('\n'.join(list_to_write))

        ''' Processing structure file if it is from PDB database '''
        if "AF-" not in structure_id:
            extract_chain(structure_path, f"{structure_id}_{chain}", preprocessing_path)
        else:
            if os.path.isfile(structure):
                shutil.copyfile(os.path.join(structure_path, basename_structure_file), os.path.join(preprocessing_path, f"{structure_id}.pdb"))
            else:
                shutil.copyfile(os.path.join(structure_path, f"{structure_id}.pdb"), os.path.join(preprocessing_path, f"{structure_id}.pdb"))

    ''' Write structures uncorrected '''
    if len(uncorrected_structure_files) > 0:
        with open(log_file, "w") as file:
            file.write('Uncorrected structure ID:\n' + '\n'.join(uncorrected_structure_files) + '\n')

''' Launch script '''
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")