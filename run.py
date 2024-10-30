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

    parser = argparse.ArgumentParser(description='### The tool for predicting probabilities of proteolytic events in proteins by structure information (structures from PDB or AlphaFoldDB databases) ###')
    parser.add_argument("-i", "--input", help="The structure ID from PDB or AlphaFoldDB databases, or a full path to the PDB file from PDB or AlphaFoldDB databases, or a full path to the file with a table where the first column is the structure ID, and the second is the name of the chain (You need specify '-' for AlphaFold structure). The value separator should be a comma. For example, 1) -i '4GAW' -c 'A' 2) -i 'AF-P00257-F1' 3) -i '/dir1/dir2/4GAW.pdb' -c 'A' 4) -i '/dir1/dir2/AF-P00257-F1.pdb' 5) -i '/dir1/dir2/structures.txt'. Also see directory 'input_examples' for a better understanding.", type=str)
    parser.add_argument("-c", "--chain", help="A name of the chain, only for PDB ID or PDB file. For example, 'A'.", type=str)
    parser.add_argument("-l", "--locally", help="Run DSSP program locally (specify '-l') or via server.", action="store_true")
    parser.add_argument("-v", "--visualise", help="Visualise the total proteolytic score (specify '-v') or not.", action="store_true")
    args = parser.parse_args()
   
    print("######## Step 1 -- Generating of PDB file with separate chain ########")
    try:
        # ID #
        if (not args.input is None) and (not os.path.isfile(args.input)):
            # AlphaFold structure #
            if 'AF-' in args.input:
                process1 = subprocess.run(["python", f"{script_path}/1.download_structure_with_sep_chain.py", f"-i={args.input}"], capture_output=True, text=True, check=True)
            # PDB structure #
            else:
                process1 = subprocess.run(["python", f"{script_path}/1.download_structure_with_sep_chain.py", f"-i={args.input}", f"-c={args.chain}"], capture_output=True, text=True, check=True)
        # PDB file #
        elif (not args.input is None) and os.path.isfile(args.input) and args.input.endswith('.pdb'):
            # AlphaFold structure #
            if 'AF-' in args.input:
                process1 = subprocess.run(["python", f"{script_path}/1.download_structure_with_sep_chain.py", f"-i={args.input}"], capture_output=True, text=True, check=True)
            # PDB structure #
            else:
                process1 = subprocess.run(["python", f"{script_path}/1.download_structure_with_sep_chain.py", f"-i={args.input}", f"-c={args.chain}"], capture_output=True, text=True, check=True)
        # TXT file #
        elif (not args.input is None) and os.path.isfile(args.input) and args.input.endswith('.txt'):
            process1 = subprocess.run(["python", f"{script_path}/1.download_structure_with_sep_chain.py", f"-i={args.input}"], capture_output=True, text=True, check=True)
        else:
            print("ERROR: Possibly, your input parameters (-i or -c) are uncorrected! Check with templates:\n1: -i 4GAW -c A\n2: -i AF-P00257-F1\n3: -i /dir1/dir2/4GAW.pdb -c A\n4: -i /dir1/dir2/AF-P00257-F1.pdb\n5: -i /dir1/dir2/structures.txt\n")
            sys.exit()
        
        if 'ERROR' in process1.stdout:
            print(process1.stdout)
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
    print("######## Step 2.1. -- Ligands and solvent removing ########")
    try:
        process2 = subprocess.run(["chimera", "--nogui", f"--script={script_path}/2.del_ligands_solvent_pyv2.py"], capture_output=True, text=True, check=True)
        if "Error while processing" in process2.stderr:
            error_index = process2.stderr.find("Error while processing")
            print(process2.stderr[error_index:])
            sys.exit()
    except subprocess.CalledProcessError as e:
        print()
        print(e.stderr)
        sys.exit()

    print("######## Step 2.2. -- PDB-file parsing ########")
    try:
        process3 = subprocess.run(["python", f"{script_path}/3.parse_pdb.py"], capture_output=True, text=True, check=True)
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
            process4 = subprocess.run(["python", f"{script_path}/4.generate_dssp_file.py", "-l"], capture_output=True, text=True, check=True)
        else:
            process4 = subprocess.run(["python", f"{script_path}/4.generate_dssp_file.py"], capture_output=True, text=True, check=True)
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
        process5 = subprocess.run(["python", f"{script_path}/5.extract_features.py"], capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print()
        print(e.stderr)
        sys.exit()
    print("######## Step 4 -- WELL DONE! ########\n")
    
    print("######## Step 5 -- Normalisation of data ########")
    try:
        process6 = subprocess.run(["python", f"{script_path}/6.norm_data.py"], capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print()
        print(e.stderr)
        sys.exit()
    print("######## Step 5 -- WELL DONE! ########\n")
    
    print("######## Step 6 -- Predicting structural scores ########")
    try:
        process7 = subprocess.run(["python", f"{script_path}/7.generate_structural_score.py"], capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print()
        print(e.stderr)
        sys.exit()
    print("######## Step 6 -- WELL DONE! ########\n")
    
    if args.visualise:
        print("######## Step 7 -- Visualising scores ########")
        try:
            process8 = subprocess.run(["chimera", "--nogui", f"--script={script_path}/8.map_scores_pyv2.py"], capture_output=True, text=True, check=True)
            if "Error while processing" in process8.stderr:
                error_index = process8.stderr.find("Error while processing")
                print(process8.stderr[error_index:])
                sys.exit()
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