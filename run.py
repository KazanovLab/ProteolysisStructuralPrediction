import shutil
import os
import subprocess
import sys
import time
import argparse

current_path = os.getcwd()
script_path = os.path.join(current_path, "scripts")
result_path = os.path.join(current_path, "results")
log_file = os.path.join(result_path, "report.log")

def main():
    print("######## START ########\n")

    parser = argparse.ArgumentParser(description='### CutProba - the tool for predicting probabilities of proteolytic events in proteins (Structural model part) ###')
    parser.add_argument("-i", "--input", help="PDB ID, or a full path to the file with a table where the first column is the PDB ID, and the second is the name of the chain. The value separator should be a comma. For example, '4GAW' or '/dir1/dir2/structures.txt'.", type=str)
    parser.add_argument("-c", "--chain", help="A name of the chain, if you specify only one PDB ID or one PDB file. For example, 'A'.", type=str)
    parser.add_argument("-f", "--file", help="A full path to the PDB file or a full path to the file with a table where the first column is a full path to the PDB file, and the second is the name of the chain. The value separator should be a comma. For example, '/dir1/dir2/4GAW.pdb' or '/dir1/dir2/structures.txt'.", type=str)
    parser.add_argument("-l", "--locally", help="Run DSSP locally (specify '-l') or via server", action="store_true")
    parser.add_argument("-v", "--visualise", help="Visualise total proteolytic score (specify '-v') or not", action="store_true")
    args = parser.parse_args()
   
    print("######## Step 1 -- Generating of PDB file with separate chain ########")
    try:
        if ((not (args.input is None)) and (not os.path.isfile(args.input))) and (not (args.chain is None)) and (args.file is None):
            process1 = subprocess.run(["python", f"{script_path}/1_download_structure_with_sep_chain.py", f"-i={args.input}", f"-c={args.chain}"], capture_output=True, text=True, check=True)
        elif ((not (args.input is None)) and os.path.isfile(args.input)) and (args.chain is None) and (args.file is None):
            process1 = subprocess.run(["python", f"{script_path}/1_download_structure_with_sep_chain.py", f"-i={args.input}"], capture_output=True, text=True, check=True)
        elif (args.input is None) and (not (args.chain is None)) and ((not (args.file is None)) and ('.pdb' in args.file)):
            process1 = subprocess.run(["python", f"{script_path}/1_download_structure_with_sep_chain.py", f"-f={args.file}", f"-c={args.chain}"], capture_output=True, text=True, check=True)
        elif (args.input is None) and (args.chain is None) and ((not (args.file is None)) and ('.pdb' not in args.file) and os.path.isfile(args.file)):
            process1 = subprocess.run(["python", f"{script_path}/1_download_structure_with_sep_chain.py", f"-f={args.file}"], capture_output=True, text=True, check=True)
        else:
            print("ERROR: Possibly, your input of PDB ID (-i), chain (-c) or PDB file (-f) is uncorrected! Check with templates:\n1: -i 4GAW -c A\n2: -i /dir1/dir2/structures.txt\n3: -f /dir1/dir2/name.pdb -c A\n4: -f /dir1/dir2/pdb_structures.txt\n")
            sys.exit()
    except subprocess.CalledProcessError as e:
        print()
        print(e.stderr)
        sys.exit()
    # Checkpoint 1 #
    if (len(os.listdir(result_path)) == 1) and (os.path.exists(log_file)):
        print("\n******** WARNING ********\nNo structures are downloaded!\n")
        sys.exit()
    print("######## Step 1 -- WELL DONE! ########\n")
    
    print("######## Step 2 -- Preprocessing  ########")
    print("######## Step 2.1. -- Ligand removing ########")
    try:
        process2 = subprocess.run(["chimera", "--nogui", f"--script={script_path}/2_del_ligands_pyv2.py"], capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print()
        print(e.stderr)
        sys.exit()
    
    print("######## Step 2.2. -- PDB-file parsing ########")
    try:
        process3 = subprocess.run(["python", f"{script_path}/3_parse_pdb.py"], capture_output=True, text=True, check=True)
        if 'ERROR' in process3.stdout:
            print(process3.stdout)
            sys.exit()
    except subprocess.CalledProcessError as e:
        print()
        print(e.stderr)
        sys.exit()
    print("######## Step 2 -- WELL DONE! ########\n")
    
    print("######## Step 3 -- Generating the DSSP file ########")
    try:
        if args.locally:
            process4 = subprocess.run(["python", f"{script_path}/4_create_dssp.py", "-l"], capture_output=True, text=True, check=True)
        else:
            process4 = subprocess.run(["python", f"{script_path}/4_create_dssp.py"], capture_output=True, text=True, check=True)
            if "ERROR" in process4.stdout:
                print(process4.stdout)
                sys.exit()
    except subprocess.CalledProcessError as e:
        print()
        print(e.stderr)
        sys.exit()
    print("######## Step 3 -- WELL DONE! ########\n")
    
    print("######## Step 4 -- Extracting features ########")
    try:
        process5 = subprocess.run(["python", f"{script_path}/5_extract_features.py"], capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print()
        print(e.stderr)
        sys.exit()
    print("######## Step 4 -- WELL DONE! ########\n")
    
    print("######## Step 5 -- Normalisation of data ########")
    try:
        process6 = subprocess.run(["python", f"{script_path}/6_norm_data.py"], capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print()
        print(e.stderr)
        sys.exit()
    print("######## Step 5 -- WELL DONE! ########\n")
    
    print("######## Step 6 -- Predicting structural scores ########")
    try:
        process7 = subprocess.run(["python", f"{script_path}/7_get_structural_score.py"], capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print()
        print(e.stderr)
        sys.exit()
    print("######## Step 6 -- WELL DONE! ########\n")
    
    if args.visualise:
        print("######## Step 7 -- Visualising scores ########")
        try:
            process8 = subprocess.run(["chimera", "--nogui", f"--script={script_path}/map_scores_pyv2.py"], capture_output=True, text=True, check=True)
            print("######## Step 7 -- WELL DONE! ########\n")
        except subprocess.CalledProcessError as e:
            print()
            print(e.stderr)
            sys.exit()
        
    print("######## FINISH ########")

if __name__ == "__main__":
    start = time.ctime()
    main()
    end = time.ctime()
    print(f"\n******** {__file__} ********\nStart: {start}\nFinish: {end}\n")