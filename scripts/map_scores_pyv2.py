import os
import time
import glob

import pandas as pd
from chimera import runCommand as run

''' Set paths '''
main_path = os.getcwd()
StructureSet_path = os.path.join(main_path, "structures")
structure_list = [i for i in os.listdir(StructureSet_path) if '.' not in i]
        
def main():
    num = 1
    for structure in structure_list:        
        structure_path = os.path.join(StructureSet_path, structure)
        score_path = os.path.join(structure_path, "scores")
        
        score_files = glob.glob(os.path.join(score_path, '*.csv'))
        if len(score_files) == 0: continue
        
        for score_file in score_files:
            structure_name = score_file.split('/')[-1].split('.csv')[0]
            print "{0} --- {1}".format(num, structure_name)
            num += 1
            
            ''' Generate attribute file '''
            score_df = pd.read_csv(score_file, dtype={"structure":str, "chain":str, "num_AA":str, "AA":str, "is_cut":int})
            structure_chain = score_df["structure"].unique()[0]
            
            score_column = [i for i in score_df.columns if 'StrModel_predict' in i][0]
            score_df[score_column + '_attribute'] = '\t:' + score_df["num_AA"] + '.' + score_df["chain"] + '\t' + score_df[score_column].map(str)
                
            attribute_file = os.path.join(score_path, structure_name + '.txt')
            save_file = os.path.join(score_path, structure_name + '.py')
                    
            with open(attribute_file, 'w') as file:
                file.write("attribute: Structural_score\nmatch mode: 1-to-1\nrecipient: residues\n" + '\n'.join(score_df[score_column + "_attribute"].tolist()))
            ''' Map score '''
            run("del")
            run("open {0}".format(os.path.join(structure_path, structure_chain + '.pdb')))
            run("defattr {}".format(attribute_file))
            run("rangecolor Structural_score 0 blue 0.5 red 1 yellow")
            run("save {}".format(save_file))
            os.remove(attribute_file)

if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print("Start: " + str(start) + '\n' + "Finish: " + str(finish) + '\n')