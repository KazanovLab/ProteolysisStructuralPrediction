""" Import libraries """
from Bio.PDB import *
import requests
import pandas as pd

import os
import sys
import glob
import time
import shutil
import argparse

''' Terminal '''
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="PDB ID, or a full path to the file with a table where the first column is the PDB ID, and the second is the name of the chain. The value separator should be a comma. For example, '4GAW' or '/dir1/dir2/structures.txt'.", type=str)
parser.add_argument("-c", "--chain", help="A name of the chain, if you specify only one PDB ID or one PDB file. For example, 'A'.", type=str)
parser.add_argument("-f", "--file", help="A full path to the PDB file or a full path to the file with a table where the first column is a full path to the PDB file, and the second is the name of the chain. The value separator should be a comma. For example, '/dir1/dir2/4GAW.pdb' or '/dir1/dir2/structures.txt'.", type=str)
args = parser.parse_args()

''' Set paths '''
main_path = os.getcwd()
StructureSet_path = os.path.join(main_path, "results")
if not os.path.exists(StructureSet_path):
    os.mkdir(StructureSet_path)

''' Logging '''    
log_file = os.path.join(StructureSet_path, "report.log")
''' List of structures '''   
if ((not (args.input is None)) and (not os.path.isfile(args.input))) and (not (args.chain is None)) and (args.file is None):
    structure_list = [[args.input, args.chain]]
elif ((not (args.input is None)) and os.path.isfile(args.input)) and (args.chain is None) and (args.file is None):
    input_data = pd.read_csv(args.input, names=["structure", "chain"], dtype={"structure":str, "chain":str})
    structure_list = input_data.values
elif (args.input is None) and (not (args.chain is None)) and ((not (args.file is None)) and ('.pdb' in args.file)):
    structure_list = [[args.file, args.chain]]
elif (args.input is None) and (args.chain is None) and ((not (args.file is None)) and ('.pdb' not in args.file) and os.path.isfile(args.file)):
    input_data = pd.read_csv(args.file, names=["structure", "chain"], dtype={"structure":str, "chain":str})
    structure_list = input_data.values
#else:
#    print("ERROR: Possibly, your input of PDB ID (-i), chain (-c) or PDB file (-f) is uncorrected! Check with templates:\n1: -i 4GAW -c A\n2: -i /dir1/dir2/structures.txt\n3: -f /dir1/dir2/name.pdb -c A\n4: -f /dir1/dir2/pdb_structures.txt\n")
#    sys.exit()
     
def download_PDB(structure_id):
	response = requests.get(f"https://files.rcsb.org/download/{structure_id}.pdb")
	return response.text

def download_AF(structure_id):
	response = requests.get(f"https://alphafold.ebi.ac.uk/files/{structure_id}-model_v1.pdb")
	return response.text

def process_pdb_file(path_to_read, structure_id, structure_chain, path_to_save):
    class ChainSelect(Select):

        def __init__(self, needed_chain):
            self.needed_chain = needed_chain

        def accept_chain(self, chain):
            if chain.id == self.needed_chain:
                return 1
            else:
                return 0

    parser = PDBParser()
    io = PDBIO()
    structure_io = parser.get_structure(structure_id, f"{path_to_read}/{structure_id}.pdb")
    if len(structure_io) > 1:
        input("CHECK!")
        sys.exit()
    else:
        io.set_structure(structure_io)
        io.save(f"{path_to_save}/{structure_id}_{structure_chain}.pdb", ChainSelect(structure_chain))

def main():
    
    uncorrected_structure_files = []
    for num, (structure, chain) in enumerate(structure_list):
        #print(f"{num + 1} --- {structure}")
        print(num, structure, chain)
        
        if os.path.isfile(structure) and ('.pdb' in structure):
            structure_file = os.path.basename(structure)
            structure_id = structure_file.split('.')[0]
        else:
            structure_id = structure
        
        structure_path = os.path.join(StructureSet_path, structure_id)
        if not os.path.exists(structure_path):
            os.mkdir(structure_path)
             
        ''' Download structure file if required'''
        if os.path.isfile(structure) and ('.pdb' in structure):
            shutil.copyfile(structure, os.path.join(structure_path, structure_file))  
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
            process_pdb_file(structure_path, structure_id, chain, structure_path)

    ''' Write structures uncorrected '''
    if len(uncorrected_structure_files) > 0:
        with open(log_file, "w") as file:
            file.write('Uncorrected structure ID\n' + '\n'.join(uncorrected_structure_files) + '\n')

''' Launch script '''
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")
