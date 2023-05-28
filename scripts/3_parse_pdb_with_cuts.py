''' Import libraries '''
import os
import sys
import glob
import time

import pandas as pd
import numpy as np

''' Set paths '''
main_path = "/data2/Protease/TEST"
dataset_path = os.path.join(main_path, "datasets")
name_of_set = input("Enter the name of structure set:\n{}\n".format(', '.join([i for i in os.listdir(dataset_path) if '.' not in i])))
substrate_path = os.path.join(main_path, "input", "pe_info", "proteolytic_events", name_of_set)
substrate_files = [i for i in os.listdir(substrate_path) if '.csv' in i]
StructureSet_path = os.path.join(dataset_path, name_of_set)
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
    path = os.path.split(structure_file)[0]
    name_structure = structure_file.split('/')[-1].split('.')[0].split('_')[0]
    name_chain = structure_file.split('/')[-1].split('.')[0].split('_')[1]
    
    with open(os.path.join(path, f"{name_structure}.pdb")) as file:
        data = file.read()
    
    seqres = [i.split() for i in data.strip('\n').split('\n') if ('SEQRES' == i.split()[0]) and (name_chain == i.split()[2])]
    len_sequence = int(seqres[0][3])
    sequence = ''.join([''.join([aminoacid_code[j] if j in aminoacid_code else 'X' for j in i[4:]]) for i in seqres])
    if len_sequence != len(sequence):
        print(sequence)
        print(len_sequence, len(sequence))
        input('CHECK!')
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

def drop_duplicate_columns(data, name_column_suffix):
    selected_columns = [i for i in data.columns if name_column_suffix in i]
    for i, selected_column in enumerate(selected_columns):
        if i == 0:
            data[name_column_suffix.lstrip('_')] = data[selected_column]
        else:
            if name_column_suffix == "_is_cut":
                data[name_column_suffix.lstrip('_')] += data[selected_column]
            else:
                data[name_column_suffix.lstrip('_')] += '&' + data[selected_column]
        del data[selected_column]
    return data
                        
def main():
    
    num = 1
    for structure in structure_list:
        
        structure_path = os.path.join(StructureSet_path, structure)
        structure_files = glob.glob(os.path.join(structure_path, "*_*.pdb")) + glob.glob(os.path.join(structure_path, "AF-*-F1.pdb"))
        
        for structure_file in structure_files:
            structure_name = structure_file.split('/')[-1].split('.')[0]
            print(f"{num} --- {structure_name}")
            num += 1
            
            ''' Get all substrates for specific structure '''
            substrate_dict = {}
            for substrate_file in substrate_files:
                substrate_name = substrate_file.split('.')[0]
                substrate_data = pd.read_csv(os.path.join(substrate_path, substrate_file), dtype={"SUBSTRATE_num_AA":str, "is_cut":int, "MEROPS_CODE":str, "SUBSTRATE_name":str, "structure":str})
                if [structure_name] == substrate_data["structure"].unique():
                    substrate_dict[substrate_name] = [''.join(substrate_data["SUBSTRATE_AA"].tolist()), substrate_data]
                if len(substrate_data["structure"].unique()) > 1:
                    print(substrate_data["structure"].unique())
                    input('CHECK!!!')
                    sys.exit()
            
            ''' Get ATOM information '''
            atom_data = pd.DataFrame(parse_atom(structure_file))
            atom_sequence = ''.join(atom_data["ATOM_AA"].tolist())
        
            if 'AF-' not in structure:
                ''' Get SEQRES information '''
                seqres_sequence, chain = get_seqres(structure_file)
            else:
                seqres_sequence, chain = atom_sequence, atom_data['chain'].unique()[0]
            
            sequence_path = os.path.join(structure_path, "sequence")
            if not os.path.exists(sequence_path):
                os.mkdir(sequence_path)
            feature_path = os.path.join(structure_path, "features")
            if not os.path.exists(feature_path):
                os.mkdir(feature_path)
            
            ''' Apply ClustalO to map SEQRES and ATOM information '''
            input_CW_file = os.path.join(sequence_path, f"{structure_name}_SUBSTRATE.fa")
            output_CW_file = os.path.join(sequence_path, f"{structure_name}_SUBSTRATE_ClustalO.fa")
            
            clustalO_input = f">SEQRES\n{seqres_sequence}\n>ATOM\n{atom_sequence}\n"
            for substrate_name in substrate_dict:
                clustalO_input += f">{substrate_name}_SUBSTRATE\n{substrate_dict[substrate_name][0]}\n"
                
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
            for substrate_name in substrate_dict:
                substrate_data = substrate_dict[substrate_name][1]
                del substrate_data["structure"]
                substrate_data["is_cut"] = substrate_data["is_cut"].astype(str)
                new_columns = {i:substrate_name + '_' + i for i in substrate_data.columns}
                substrate_data = substrate_data.rename(columns=new_columns)

                SUBSTRATE_AA_index = map_data.loc[map_data[f"{substrate_name}_SUBSTRATE_AA"] != '-', f"{substrate_name}_SUBSTRATE_AA"].index
                substrate_data = substrate_data.set_index(SUBSTRATE_AA_index)
                
                final_data = pd.concat([final_data, substrate_data[[i for i in substrate_data.columns if i != f"{substrate_name}_SUBSTRATE_AA"]]], axis=1)
            
            ''' Preprocessing final data '''
            final_data = final_data.fillna('-')
            final_data["structure"] = structure_name
            final_data = drop_duplicate_columns(final_data, '_MEROPS_CODE')
            final_data = drop_duplicate_columns(final_data, '_is_cut')
            final_data = drop_duplicate_columns(final_data, '_SUBSTRATE_name')
            final_data["is_cut"] = final_data["is_cut"].apply(lambda x: 1 if ('1' in x) else 0) 
            
            final_columns = ["structure", "chain", "ATOM_AA", "ATOM_num_AA", "bfac", "is_cut", "MEROPS_CODE", "SUBSTRATE_name"]
            final_atom_data = final_data.loc[final_data["ATOM_AA"] != '-', final_columns]
            if 'AF-' in structure_name:
                final_atom_data = final_atom_data.rename(columns={"bfac":"AF_score"})
            final_atom_data = final_atom_data.rename(columns={"ATOM_AA":"AA", "ATOM_num_AA":"num_AA"})
            
            ''' Save results '''
            final_data.to_csv(os.path.join(feature_path, f"{structure_name}_init.csv"), index=False)
            final_atom_data.to_csv(os.path.join(feature_path, f"{structure_name}.csv"), index=False)
            
''' Launch script '''                
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")