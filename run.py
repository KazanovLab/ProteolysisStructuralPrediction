import shutil
import os
import time
import argparse

current_path = os.getcwd()
script_path = os.path.join(current_path, "scripts")

def main():
    print("######## START ########")

    parser = argparse.ArgumentParser(description='### CutProba - the tool for predicting probabilities of proteolytic events in proteins (Structural model part) ###')
    parser.add_argument("-i", "--input", help="The name of your PDB ID with separate chain, or TXT file with such structures separated by 'ENTER', for example '4GAW_A'", type=str)
    parser.add_argument("-l", "--locally", help="Run DSSP locally (specify '-l') or via server", action="store_true")
    parser.add_argument("-v", "--visualise", help="Visualise total proteolytic score (specify '-v') or not", action="store_true")
    args = parser.parse_args()

    PDB_ID = args.input
    if '.txt' in PDB_ID:
        PDB_ID = os.path.join(current_path, PDB_ID)
    
    print("######## Step 1 -- Creating of PDB file with separate chain ########")
    os.system(f"python {script_path}/1_download_structure_with_sep_chain.py -i {PDB_ID}")
    
    print("######## Step 2 -- Preprocessing  ########")
    os.system(f"chimera --nogui --script {script_path}/2_del_ligands_pyv2.py")
    os.system(f"python {script_path}/3_parse_pdb.py")
    
    print("######## Step 3 -- Creating the DSSP file ########")
    if args.locally:
        os.system(f"python {script_path}/4_create_dssp.py -l")
    else:
        os.system(f"python {script_path}/4_create_dssp.py")
    
    print("######## Step 4 -- Extracting features ########")
    os.system(f"python {script_path}/5_extract_features.py")
    
    print("######## Step 5 -- Normalisation of data ########")
    os.system(f"python {script_path}/6_norm_data.py")
    
    print("######## Step 6 -- Predicting structural scores ########")
    os.system(f"python {script_path}/7_get_structural_score.py")
    
    if args.visualise:
        print("######## Step 7 -- Visualising scores ########")
        os.system(f"chimera --nogui --script {script_path}/map_scores_pyv2.py")

    print("######## FINISH ########")

if __name__ == "__main__":
    start = time.ctime()
    main()
    end = time.ctime()
    print(f"\n******** {__file__} ********\nStart: {start}\nFinish: {end}\n")