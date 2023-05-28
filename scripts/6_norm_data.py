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

''' Logging '''
log_file = os.path.join(StructureSet_path, "report.log")

def normalise(dataset, name_input_vector, name_output_vector):
    min_value = np.min(dataset[name_input_vector])
    max_value = np.max(dataset[name_input_vector])
    if max_value - min_value == 0:
        dataset[name_output_vector] = np.nan
    else:
        dataset[name_output_vector] = dataset[name_input_vector].apply(lambda x: round( (x - min_value) / (max_value - min_value), 3))

def main():
    warnings = []
    
    num = 1
    for structure in structure_list:        
        structure_path = os.path.join(StructureSet_path, structure)
        feature_path = os.path.join(structure_path, "features")
        
        feature_files = glob.glob(os.path.join(feature_path, '*_final.csv'))
        for feature_file in feature_files:
            structure_name = '_'.join(feature_file.split('/')[-1].split('.')[0].split('_')[:-1])
            #print(f"{num} --- {structure_name}")
            num += 1
            
            feature_df = pd.read_csv(feature_file)
            if ('AF-' in structure_name) and ('is_cut' in feature_df.columns):
                feature_df = feature_df.astype({"structure":str, "chain":str, "AA":str, "num_AA":str, "is_cut":int, "MEROPS_CODE":str, "SUBSTRATE_name":str, "init_SS_type":str, "SS_type":str, "bp1":int, "bp2":int, "ACC":float, "AF_score":float, "len_loop":int, "is_terminus":int})
            elif ('AF-' not in structure_name) and ('is_cut' in feature_df.columns):
                feature_df = feature_df.astype({"structure":str, "chain":str, "AA":str, "num_AA":str, "is_cut":int, "MEROPS_CODE":str, "SUBSTRATE_name":str, "init_SS_type":str, "SS_type":str, "bp1":int, "bp2":int, "ACC":float, "bfac":float, "len_loop":int, "is_terminus":int})
            elif ('AF-' in structure_name) and ('is_cut' not in feature_df.columns):
                feature_df = feature_df.astype({"structure":str, "chain":str, "AA":str, "num_AA":str, "init_SS_type":str, "SS_type":str, "bp1":int, "bp2":int, "ACC":float, "AF_score":float, "len_loop":int, "is_terminus":int})
            else:
                feature_df = feature_df.astype({"structure":str, "chain":str, "AA":str, "num_AA":str, "init_SS_type":str, "SS_type":str, "bp1":int, "bp2":int, "ACC":float, "bfac":float, "len_loop":int, "is_terminus":int})
            
            ''' ACC (structure) '''
            normalise(feature_df, "ACC", "ACC_N_chain")
            if np.unique(feature_df["ACC_N_chain"].isna()) == [True]:
                warnings.append(f"{structure_name} - uncorrected normalisation of Solvent accessibility (ACC is constant)!\n******** WARNING ********\nPredictions may be uncorrected!")
                feature_df["ACC_N_chain"].fillna(0, inplace=True)
                #continue
            ''' B-factor / AlphaFold prediction score '''
            if 'AF-' in structure_name:
                normalise(feature_df, "AF_score", "AF_score_N_chain")
                if np.unique(feature_df["AF_score_N_chain"].isna()) == [True]:
                    warnings.append(f"{structure_name} - uncorrected normalisation of AlphaFold prediction score because of values are same!\n******** WARNING ********\nPredictions may be uncorrected!")
                    feature_df["AF_score_N_chain"].fillna(0, inplace=True)
                    #continue
            else:
                normalise(feature_df, "bfac", "bfac_N_chain")
                if np.unique(feature_df["bfac_N_chain"].isna()) == [True]:
                    warnings.append(f"{structure_name} - uncorrected normalisation of B-factor! B-factor is constant!\n******** WARNING ********\nPredictions may be uncorrected!")
                    feature_df["bfac_N_chain"].fillna(0, inplace=True)
                    #continue
            ''' Loop length '''
            normalise(feature_df, "len_loop", "len_loop_N_chain")
            if np.unique(feature_df["len_loop_N_chain"].isna()) == [True]:
                warnings.append(f"{structure_name} - uncorrected normalisation of Loop length because of values are same. Structure has only one type of secondary structure!\n******** WARNING ********\nPredictions may be uncorrected!")
                feature_df["len_loop_N_chain"].fillna(0, inplace=True)
                #continue

            ''' Secondary structure type coding '''
            feature_df = pd.concat([pd.get_dummies(feature_df, prefix='SS', columns=['SS_type']), feature_df["SS_type"]], axis=1)
            if 'SS_O' not in feature_df.columns:
                feature_df["SS_O"] = 0
            if 'SS_B' not in feature_df.columns:
                feature_df["SS_B"] = 0
            if 'SS_H' not in feature_df.columns:
                feature_df["SS_H"] = 0
            
            ''' Save results '''
            feature_df.to_csv(os.path.join(feature_path, feature_file.split('.csv')[0] + '_norm.csv'), index=False)
    
    if len(warnings) > 0:
        if os.path.exists(log_file):
            with open(log_file, 'a') as file:
                file.write('Uncorrected normalisation\n' + '\n'.join(warnings) + '\n')
        else:
            with open(log_file, 'w') as file:
                file.write('Uncorrected normalisation\n' + '\n'.join(warnings) + '\n')
                
''' Launch script '''                
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")