''' Import libraries '''
import pandas as pd
import numpy as np

import os
import sys
import glob
import time

''' Set paths '''
main_path = os.getcwd()
StructureSet_path = os.path.join(main_path, "results")
structure_list = [i for i in os.listdir(StructureSet_path) if '.' not in i]

''' Set AA code '''
aminoacid_code = {'GLY':'G',
                  'ALA':'A',
                  'VAL':'V',
                  'LEU':'L',
                  'ILE':'I',
                  'CYS':'C',
                  'MET':'M',
                  'PHE':'F',
                  'TYR':'Y',
                  'TRP':'W',
                  'PRO':'P',
                  'SER':'S',
                  'THR':'T',
                  'ASN':'N',
                  'GLN':'Q',
                  'ASP':'D',
                  'GLU':'E',
                  'HIS':'H',
                  'LYS':'K',
                  'ARG':'R'
                  }

def get_seqres(structure_file):
    path = os.path.split(os.path.split(structure_file)[0])[0]
    name_structure = '_'.join(structure_file.split('/')[-1].split('.del.pdb')[0].split('_')[:-1])
    name_chain = structure_file.split('/')[-1].split('.del.pdb')[0].split('_')[-1]
    
    with open(os.path.join(path, f"{name_structure}.pdb")) as file:
        data = file.read()
    
    seqres = [i.split() for i in data.strip('\n').split('\n') if ('SEQRES' == i.split()[0]) and (name_chain == i.split()[2])]
    if len(seqres) == 0:
        print(f"ERROR: Check chain in {structure_file}, please! Possibly, this chain doesn't exist!")
        sys.exit()
    len_sequence = int(seqres[0][3])
    sequence = ''.join([''.join([aminoacid_code[j] if j in aminoacid_code else 'X' for j in i[4:]]) for i in seqres])
    if len_sequence != len(sequence):
        print(sequence)
        print(len_sequence, len(sequence))
        print('CHECK! (for developers)')
        sys.exit()
        
    return sequence, name_chain

def parse_atom(structure_file):
    with open(structure_file, 'r') as file:
        data = file.read()
    
    dict_data = {}
    
    for stroka in data.strip('\n').split('\n'):
        if stroka[:6] in ["ATOM  ", "HETATM"]:
            residue = stroka[17:20]
            residue = aminoacid_code[residue] if residue in aminoacid_code else 'X' # X - modified residues
            chain = stroka[21]
            num_AA = stroka[22:27].strip() # + insertion code
            bfac = float(stroka[60:66].strip())
        
            info = f"{residue}_{chain}_{num_AA}"
            if info not in dict_data:
                dict_data[info] = [bfac]
            else:
                dict_data[info].append(bfac)
    
    list_data = [i.split('_') + [round(np.mean(dict_data[i]), 2)] for i in dict_data]
    
    dict_data = {}
    dict_data["ATOM_AA"] = []
    dict_data["chain"] = []
    dict_data["ATOM_num_AA"] = []
    dict_data["bfac"] = []
    
    for e in list_data:
        dict_data["ATOM_AA"].append(e[0])
        dict_data["chain"].append(e[1])
        dict_data["ATOM_num_AA"].append(e[2])
        dict_data["bfac"].append(e[3])
    
    return dict_data

def map_sequence(align_file):
    
    with open(align_file) as file:
        data = file.read()
    
    data = data.lstrip('>').rstrip('\n').split('>')
    map_results = {i.split('\n')[0] + '_AA':list(''.join(i.split('\n')[1:])) for i in data}
    
    map_data = pd.DataFrame(map_results)
    return map_data

def main():
    
    num = 1
    for structure in structure_list:
        
        structure_path = os.path.join(StructureSet_path, structure)
        preprocessing_path = os.path.join(structure_path, "preprocessing")
        structure_files = glob.glob(os.path.join(preprocessing_path, "*.del.pdb"))
        
        for structure_file in structure_files:
            structure_name = structure_file.split('/')[-1].split('.del.pdb')[0]
            #print(f"{num} --- {structure_name}")
            num += 1
        
            ''' Get ATOM information '''
            atom_data = pd.DataFrame(parse_atom(structure_file))
            atom_sequence = ''.join(atom_data["ATOM_AA"].tolist())
        
            if 'AF-' not in structure_name:
                ''' Get SEQRES information '''
                seqres_sequence, chain = get_seqres(structure_file)
            else:
                seqres_sequence, chain = atom_sequence, atom_data['chain'].unique()[0]
            
            sequence_path = os.path.join(structure_path, "alignments")
            if not os.path.exists(sequence_path):
                os.mkdir(sequence_path)
            feature_path = os.path.join(structure_path, "features")
            if not os.path.exists(feature_path):
                os.mkdir(feature_path)
            
            ''' Apply ClustalO to map SEQRES and ATOM information '''
            input_CW_file = os.path.join(sequence_path, f"{structure_name}.SEQRES_ATOM.fa")
            output_CW_file = os.path.join(sequence_path, f"{structure_name}.SEQRES_ATOM_ClustalO.fa")
            
            clustalO_input = f">SEQRES\n{seqres_sequence}\n>ATOM\n{atom_sequence}\n"
            
            with open(input_CW_file, 'w') as file:
                file.write(clustalO_input)
            
            os.system(f'''clustalo \
                      --infile={input_CW_file} \
                      --infmt=fa \
                      --outfile={output_CW_file} \
                      --outfmt=fa
                      ''')
                
            ''' Mapping '''
            map_data = map_sequence(output_CW_file)
            
            ATOM_AA_index = map_data.loc[map_data["ATOM_AA"] != '-', 'ATOM_AA'].index
            atom_data = atom_data.set_index(ATOM_AA_index)
            final_data = pd.concat([map_data, atom_data[[i for i in atom_data.columns if i != "ATOM_AA"]]], axis=1)
            
            final_data = final_data.fillna('-')
            final_data["structure"] = structure_name
            final_columns = ["structure", "chain", "ATOM_AA", "ATOM_num_AA", "bfac"]
            final_atom_data = final_data.loc[final_data["ATOM_AA"] != '-', final_columns]
            if 'AF-' in structure_name:
                final_atom_data = final_atom_data.rename(columns={"bfac":"confidence_score"})
            final_atom_data = final_atom_data.rename(columns={"ATOM_AA":"AA", "ATOM_num_AA":"num_AA"})
          
            ''' Save results '''
            final_data.to_csv(os.path.join(feature_path, f"{structure_name}.init.csv"), index=False)
            final_atom_data.to_csv(os.path.join(feature_path, f"{structure_name}.csv"), index=False)

''' Launch script '''                
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")