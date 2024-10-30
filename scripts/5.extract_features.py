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

def parse_dssp(dssp_file):
    SS = ['a','b','c','d','e','f','g','h',
	'i','j','k','l','m','n','o','p','q',
	'r','s','t','u','v','w','z','x','y']
    HTR = 'X'
    b_sheet = ['E']
    a_helix = ['H', 'G', 'I']
    loop = ['T', 'S', 'B']
    
    with open(dssp_file) as file:
        data = file.read()
    
    ''' DSSP columns index '''
    for num_string, stroka in enumerate(data.strip('\n').split('\n')):
        if stroka.split()[0] == '#':
            break
    
    dict_data = {}
    dict_data['num_AA'] = []
    dict_data['chain'] = []
    dict_data['AA'] = []
    dict_data['init_SS_type'] = []
    dict_data['SS_type'] = []
    dict_data['bp1'] = []
    dict_data['bp2'] = []
    dict_data['ACC'] = []
    
    for stroka in data.strip('\n').split('\n')[(num_string + 1):]:
        num_AA = stroka[5:11].strip()
        chain = stroka[11].strip()
        residue = stroka[12:14].strip()
        #residue = 'C' if residue in SS else residue
        init_SS_type = 'U' if stroka[15:17].strip() == '' else stroka[15:17].strip()
        SS_type = 'B' if init_SS_type in b_sheet else 'H' if init_SS_type in a_helix else 'O'
        bp1 = 1 if int(stroka[25:29].strip()) > 0 else 0
        bp2 = 1 if int(stroka[29:33].strip()) > 0 else 0 # 33 - letter code of beta sheet
        ACC = float(stroka[34:38].strip())
        
        if residue == '!':
            continue
        elif residue in SS:
            print(stroka)
            print(ord_num, num_AA, chain, residue, init_SS_type, SS_type)
            print('1 - CHECK! (for developers)')
            sys.exit()
        if len(residue) > 1:
            print(stroka)
            print(ord_num, num_AA, chain, residue, init_SS_type, SS_type)
            print('2 - CHECK! (for developers)')
            sys.exit()
        if len(init_SS_type) > 1:
            print(stroka)
            print(ord_num, num_AA, chain, residue, init_SS_type, SS_type)
            print('3 - CHECK! (for developers)')
            sys.exit()
        
        dict_data['num_AA'].append(num_AA)
        dict_data['chain'].append(chain)
        dict_data['AA'].append(residue)
        dict_data['init_SS_type'].append(init_SS_type)
        dict_data['SS_type'].append(SS_type)
        dict_data['bp1'].append(bp1)
        dict_data['bp2'].append(bp2)
        dict_data['ACC'].append(ACC)
    
    return dict_data

def get_len_loop(SS_type):
    
    len_loop = []
    indices_O = []
    count_O = 0
    for i, SS in enumerate(SS_type):
        len_loop.append(0)
        if SS == 'O':
            count_O += 1
            indices_O.append(i)
        else:   
            for j in indices_O:
                len_loop[j] = count_O
            count_O = 0
            indices_O = []
    for j in indices_O: # for C-terminus loop
        len_loop[j] = count_O
    
    return len_loop

def get_termini(init_SS_type):
    ''' part 1 - Structured and unstructured regions'''
    unstructured = ["U", "S", "T", "B"]
    structured = ["E", "H", "G", "I"]
    parts_STR = []
    local_part = ''
    index = 0

    while index != len(init_SS_type) - 1:
        value = init_SS_type[index]
        next_value = init_SS_type[index + 1]
        if (value in unstructured and next_value in unstructured) or (value in structured and next_value in structured):
            local_part += value
        else:
            local_part += value
            parts_STR.append(local_part)
            local_part = ''
        index += 1
    if local_part == '':
        parts_STR.append(init_SS_type[-1])
    else:
        parts_STR.append(local_part + init_SS_type[-1])
    
    ''' part 2 - '''
    def process_STR(parts_STR):
        bends = ["SS", "TT", "ST", "TS", "US", "UT", "SU", "TU", "BS", "SB", "BT", "TB"]
        edge = ''

        # part 1
        part_1 = parts_STR[0]
        #print(f"Part 1: {part_1}")
        bend_counts = 0

        if len(part_1) % 2 == 0:
            for index in range(0, len(part_1), 2):
                edge += part_1[index] + part_1[index + 1]
                if part_1[index] + part_1[index + 1] in bends:
                    bend_counts += 1
                if bend_counts == 3:
                    #print(f"\nBend counts: {bend_counts}")
                    #print(f"Edge: {edge}")
                    return edge
        else:
            if len(part_1) == 1:
                edge = part_1
            else:
                for index in range(len(part_1) - 1):
                    if index % 2 == 0:
                        edge += part_1[index] + part_1[index + 1]
                        if part_1[index] + part_1[index + 1] in bends:
                            bend_counts += 1
                    elif index == len(part_1) - 2:
                        edge += part_1[index + 1]
                        if part_1[index] + part_1[index + 1] in bends:
                            bend_counts += 1

                    if bend_counts == 3:
                        #print(f"\nBend counts: {bend_counts}")
                        #print(f"Edge: {edge}")
                        return edge

        # part 2
        part_2 = parts_STR[1]
        part_3 = parts_STR[2]
        #print(f"Part 2: {part_2}")
        #print(f"Part 3: {part_3}")
        
        if len(set(part_2)) > 1:
            edge += part_2
        else:
            # part 3
            if len(part_3) < 3:
                # Split part 4 or not? #
                if (part_3 in bends) or (len(parts_STR) == 3):
                    edge += part_2 + part_3
                else:
                    part_4 = parts_STR[3]
                    #print(f"Part 4: {part_4}")
                    if len(set(part_4)) > 1:
                        if ("H" in set(part_4) and "I" in set(part_4)) or ("H" in set(part_4) and "G" in set(part_4)) or ("G" in set(part_4) and "I" in set(part_4)):
                            edge += part_2 + part_3 + part_4
                        else:
                            for i, l in enumerate(part_4):
                                if l != part_4[i + 1]:
                                    edge += part_2 + part_3 + part_4[:(i + 1)]
                                    break
                    else:
                        edge += part_2 + part_3 + part_4
            else:
                if part_3[:2] in bends:
                    edge += part_2 + part_3[:2]
                elif part_3[:3] in ["UUU", "UUS", "UUT", "UUB", "BUS", "BUT", "BUU", "UBS", "UBT", "UBU", "UBB", "BUB", "BBU", "BBB", "BBS"]:
                    edge += part_2 + part_3[:3]
                else:
                    print(part_3[:3])
                    print("4 - CHECK! (for developers)")
                    sys.exit()
                    
        #print(f"Edge: {edge}")
        return edge
    
    if len(parts_STR) < 3:
        #print(parts_STR)
        #input('5 - CHECK! (for developers)')
        termini_feature = [-1] * len(init_SS_type)
    else:
        N_terminus = process_STR(parts_STR)
        C_terminus = process_STR([i[::-1] for i in parts_STR[::-1]])
        if len(init_SS_type) - len(N_terminus) - len(C_terminus) < 0:
            #input('6 - CHECK! (for developers)')
            termini_feature = [1] * len(init_SS_type)
        else:
            termini_feature = [1] * len(N_terminus) + [0] * (len(init_SS_type) - len(N_terminus) - len(C_terminus)) + [1] * len(C_terminus)
    
    return termini_feature

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
        feature_path = os.path.join(structure_path, "features")
        sequence_path = os.path.join(structure_path, "sequence")
        
        dssp_files = glob.glob(os.path.join(structure_path, '*.dssp'))
        for dssp_file in dssp_files:
            structure_name = dssp_file.split('/')[-1].split('.')[0]
            #print(f"{num} --- {structure_name}")
            num += 1
            
            ''' ACC, bp1, bp2 '''
            dssp_data = pd.DataFrame(parse_dssp(dssp_file))
            dssp_data.loc[dssp_data['ACC'] > 200, 'ACC'] = 200
            
            ''' Length of loops feature '''
            dssp_data['len_loop'] = get_len_loop(dssp_data['SS_type'].tolist())
            dssp_data.loc[dssp_data['len_loop'] > 20, 'len_loop'] = 20
            
            ''' Termini '''
            dssp_data["is_terminus"] = get_termini(dssp_data['init_SS_type'].tolist())
            
            ''' Apply ClustalO to map PDB and DSSP information '''
            pdb_data = pd.read_csv(os.path.join(feature_path, f"{structure_name}.csv"), dtype={'chain':str, 'num_AA':str})
            
            pdb_sequence = ''.join(pdb_data["AA"].tolist())
            dssp_sequence = ''.join(dssp_data["AA"].tolist())
            
            input_CW_file = os.path.join(sequence_path, f"{structure_name}.PDB_DSSP.fa")
            output_CW_file = os.path.join(sequence_path, f"{structure_name}.PDB_DSSP_ClustalO.fa")
            
            clustalO_input = f">PDB\n{pdb_sequence}\n>DSSP\n{dssp_sequence}\n"
                
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
            
            DSSP_AA_index = map_data.loc[map_data["DSSP_AA"] != '-', 'DSSP_AA'].index
            dssp_data = dssp_data.set_index(DSSP_AA_index)
            final_data = pd.concat([map_data, pdb_data[[i for i in pdb_data.columns if i not in ["AA", "num_AA", "chain"]]], dssp_data[[i for i in dssp_data.columns if i not in ["AA"]]]], axis=1)
            
            final_data = final_data.fillna('-')
            del final_data["PDB_AA"]
            final_data = final_data.rename(columns={"DSSP_AA":"AA"})
            
            if 'AF-' in structure_name:
                final_columns = ["structure", "chain", "AA", "num_AA", "init_SS_type", "SS_type", "bp1", "bp2", "ACC", "confidence_score", "len_loop", "is_terminus"]
            else:
                final_columns = ["structure", "chain", "AA", "num_AA", "init_SS_type", "SS_type", "bp1", "bp2", "ACC", "bfac", "len_loop", "is_terminus"]
            final_data = final_data.loc[final_data["AA"] != '-', final_columns].astype({"chain":str, "num_AA":str, "bp1":int, "bp2":int, "len_loop":int, "is_terminus":int})
            
            final_data.to_csv(os.path.join(feature_path, f"{structure_name}.final.csv"), index=False)
            
''' Launch script '''                
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")