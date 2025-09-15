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

    parser = argparse.ArgumentParser(description='### The tool for predicting probabilities of proteolytic events in proteins by structure information (structures are extracted from PDB or AlphaFoldDB databases) ###')
    parser.add_argument("-s", "--structure", help="The four-symbol structure ID from PDB database (case-insensitive), or a full path to the structure file with '.pdb' extension. It must be used in conjunction with the '-c' parameter. For example, 1) -s '4GAW' -c 'A' 2) -s '/dir1/my_experimental_structure.pdb' -c 'A'. Also see directory 'input_examples' for a better understanding.", type=str)
    parser.add_argument("-m", "--model", help="The UniProt ID of a protein for downloading model structure from AlphaFoldDB database, or a full path to the model structure file with '.pdb' extension. For example, 1) -m 'P00257' 2) -m '/dir1/dir2/my_model_structure.pdb'. Also see directory 'input_examples' for a better understanding.", type=str)
    parser.add_argument("-b", "--batch", help="A batch option to generate predictions for a set of structures. A full path to the TXT-file with a table where the first column is the structure ID (The name of your AlphaFold structure file must be necessarily started with 'AF-' and the name of your AlphaFold structure ID must match a template as 'AF-{UniProt_ID}_F1'), and the second is the name of the protein chain (You need specify '-' for AlphaFold structure). The value separator should be a comma. For example, 1) -b '/dir1/dir2/structures.txt'. Also see directory 'input_examples' for a better understanding.", type=str)
    parser.add_argument("-c", "--chain", help="A name of the protein chain (case-sensitive), only used with '-s' parameter for PDB ID or PDB structure file. For example, 'A'.", type=str)
    parser.add_argument("-l", "--locally", help="Run DSSP program locally (specify '-l') or via server.", action="store_true")
    parser.add_argument("-v", "--visualise", help="Visualise the total proteolytic score (specify '-v') or not.", action="store_true")
    args = parser.parse_args()

    print("######## Step 1 -- Generating of PDB file with separate chain ########")
    try:
        
        # PDB #
        if (not args.structure is None) and (args.model is None) and (args.batch is None):
            # ID #
            if not os.path.isfile(args.structure):
                # -s 4gaw.pdb
                if args.structure.endswith('.pdb'):
                    print(f"ERROR: Your file '{args.structure}' is not found.")
                    sys.exit()
                # -s abracadabra
                if len(args.structure) != 4:
                    print(f"ERROR: The structure ID from PDB database should be four-symbol. Please, see examples in PDB database.")
                    sys.exit()
                # -s 4GAW
                if args.chain is None:
                    print(f"ERROR: Possibly, you forgot specify the specific protein chain of the structure to be analyzed. Please, specify the chain. For example, '-c A'.")
                    sys.exit()
                # -s 4GAW -c A
                process1 = subprocess.run(["python", f"{script_path}/1.download_structure_with_sep_chain.py", f"-s={args.structure}", f"-c={args.chain}"], capture_output=True, text=True, check=True)
            # File #
            else:
                # -s 4gaw.cif -c A
                if not args.structure.endswith('.pdb'):
                    print(f"ERROR: Your file '{args.structure}' doesn't have '.pdb' extension.")
                    sys.exit()
                # -s AF-124.pdb -c A
                if 'AF-' in args.structure:
                    print(f"ERROR: Your experimental structure file '{args.structure}' contains 'AF-' substring. It might be confused for generating predictions. Please, rename your file and try again!")
                    sys.exit()
                # -s /dir1/my_experimental_structure.pdb
                if args.chain is None:
                    print(f"ERROR: Possibly, you forgot specify the specific protein chain of the structure to be analyzed. Please, specify the chain. For example, '-c A'.")
                    sys.exit()
                # -s /dir1/my_experimental_structure.pdb -c A
                process1 = subprocess.run(["python", f"{script_path}/1.download_structure_with_sep_chain.py", f"-s={args.structure}", f"-c={args.chain}"], capture_output=True, text=True, check=True)
                
        # AlphaFold #
        elif (args.structure is None) and (not args.model is None) and (args.batch is None):
            # ID #
            if not os.path.isfile(args.model):
                # -m my_model_protein.pdb
                if args.model.endswith('.pdb'):
                    print(f"ERROR: Your file '{args.model}' is not found.")
                    sys.exit()
                # -m P00027
                if args.chain is None:
                    process1 = subprocess.run(["python", f"{script_path}/1.download_structure_with_sep_chain.py", f"-m={args.model}"], capture_output=True, text=True, check=True)
                # -m P00027 -c A
                else:
                    print(f"WARNING: it's not necessarily to specify '-c' parameter for model structure. The value of this parameter is not taken into account.")
                    process1 = subprocess.run(["python", f"{script_path}/1.download_structure_with_sep_chain.py", f"-m={args.model}"], capture_output=True, text=True, check=True)
            # File #
            else:
                # -m /dir1/dir2/my_model_protein.cif
                if not args.model.endswith('.pdb'):
                    print(f"ERROR: Your file '{args.model}' doesn't have 'PDB' extension.")
                    sys.exit()
                # -m /dir1/dir2/my_model_protein.pdb 
                if args.chain is None:
                    process1 = subprocess.run(["python", f"{script_path}/1.download_structure_with_sep_chain.py", f"-m={args.model}"], capture_output=True, text=True, check=True)
                # -m /dir1/dir2/my_model_protein.pdb -c A
                else:
                    print(f"WARNING: it's not necessarily to specify '-c' parameter for model structure. The value of this parameter is not taken into account.")
                    process1 = subprocess.run(["python", f"{script_path}/1.download_structure_with_sep_chain.py", f"-m={args.model}"], capture_output=True, text=True, check=True)
        
        # Batch option #
        elif (args.structure is None) and (args.model is None) and (not args.batch is None):
            # /dir1/dir2/mystructures.txt - but doesn't exist
            if not os.path.isfile(args.batch):
                print(f"ERROR: Your batch file '{args.batch}' is not found.")
                sys.exit()
            # /dir1/dir2/mystructures.tsv
            if not args.batch.endswith('.txt'):
                print(f"ERROR: Your batch file '{args.batch}' doesn't have '.txt' extension.")
                sys.exit()
            
            process1 = subprocess.run(["python", f"{script_path}/1.download_structure_with_sep_chain.py", f"-b={args.batch}"], capture_output=True, text=True, check=True)
        
        # Another possible variants #
        else:
            print("ERROR: Possibly, your input parameters are uncorrected! Check with next templates:\n1: -s 4GAW -c A\n2: -s /dir1/my_experimental_structure.pdb -c A\n3: -m P00257\n4: -m /dir1/dir2/my_model_structure.pdb\n5: -b /dir1/dir2/my_structures.txt\n")
            sys.exit()
        
        # Check stdout of the first script #
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
        print(process7.stdout)
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
