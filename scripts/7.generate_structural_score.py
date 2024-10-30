''' Import libraries '''
import pandas as pd
import numpy as np
np.random.seed(8)

import _pickle as cPickle
import os
import glob
import time

''' Set paths ''' 
main_path = os.getcwd()
model_path = os.path.join(main_path, "models")
StructureSet_path = os.path.join(main_path, "results")
structure_list = [i for i in os.listdir(StructureSet_path) if '.' not in i]

''' Features '''
features_exp = {"ACC_N_chain":float, "bfac_N_chain":float, "SS_B":int, "SS_H":int, "SS_O":int, "len_loop_N_chain":float, "is_terminus":int}
features_AF = {"ACC_N_chain":float, "confidence_score_N_chain":float, "SS_B":int, "SS_H":int, "SS_O":int, "len_loop_N_chain":float, "is_terminus":int}

''' Loading  Models '''
# Structural model #
with open(os.path.join(model_path, "structural_model.pkl"), 'rb') as file:
    structural_model = cPickle.load(file)
# Structural model for AlphaFold structures #
with open(os.path.join(model_path, "structural_model_AlphaFold.pkl"), 'rb') as file:
    structural_model_AF = cPickle.load(file)

def main():
    
    num = 1
    for structure in structure_list:        
        structure_path = os.path.join(StructureSet_path, structure)
        feature_path = os.path.join(structure_path, "features")
        
        feature_files = glob.glob(os.path.join(feature_path, '*.norm.csv'))
        for feature_file in feature_files:
            structure_name = feature_file.split('/')[-1].split('.')[0]
            #print(f"{num} --- {structure_name}")
            num += 1
            
            ''' Structural score '''
            if 'AF-' in structure_name:
                feature_df = pd.read_csv(feature_file, dtype=features_AF)
                X_test = feature_df[features_AF.keys()].values
                structural_scores = structural_model_AF.predict_proba(X_test)[:, 1]
            else:
                feature_df = pd.read_csv(feature_file, dtype=features_exp)
                X_test = feature_df[features_exp.keys()].values
                structural_scores = structural_model.predict_proba(X_test)[:, 1]

            score_path = os.path.join(structure_path, "scores")
            if not os.path.exists(score_path):
                os.mkdir(score_path)
            
            structural_score_col = f"PredictedProbability"
            feature_df[structural_score_col] = structural_scores
            feature_df[structural_score_col] = feature_df[structural_score_col].map(lambda x: round(x, 4))
            
            feature_df.to_csv(os.path.join(score_path, f"{structure_name}.StrProba.csv"), index=False)
            feature_df[["structure", "chain", "AA", "num_AA", "PredictedProbability"]].to_csv(os.path.join(structure_path, f"{structure_name}.predictions.csv"), index=False)

''' Launch script '''                
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")