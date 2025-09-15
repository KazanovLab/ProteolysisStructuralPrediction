''' Import libraries '''
import os
import sys
import glob
import time
import requests
import json
import argparse

''' Terminal '''
parser = argparse.ArgumentParser()
parser.add_argument("-l", "--locally", help="Run DSSP program locally (specify '-l') or via server", action="store_true")
args = parser.parse_args()

''' Set paths '''
main_path = os.getcwd()
StructureSet_path = os.path.join(main_path, "results")
structure_list = [i for i in os.listdir(StructureSet_path) if '.' not in i]

''' Logging '''
log_file = os.path.join(StructureSet_path, "report.log")

def download_dssp(structure_file):
    dssp_url = 'https://www3.cmbi.umcn.nl/xssp/'
    url_create = '{}api/create/pdb_file/dssp/'.format(dssp_url)
    
    try:
        response_create = requests.post(url_create, files={'file_':open(structure_file, 'rb')}, timeout=30)
    except requests.Timeout:
        print("ERROR: Timeout for response from DSSP server! Please, install DSSP and repeat this step locally!")
        sys.exit()
    except requests.ConnectionError:
        print("ERROR: No connection with DSSP server! Please, repeat later!")
        sys.exit()
        
    response_create.raise_for_status

    job_id = json.loads(response_create.text)['id']
    print("Job submitted successfully. Id is: '{}'".format(job_id))

    ready = False
    while not ready:
        local_time = time.time()
        print(int(local_time))
        url_status = '{}/api/status/pdb_file/dssp/{}/'.format(dssp_url, job_id)
        response_status = requests.get(url_status)
        response_status.raise_for_status

        status = json.loads(response_status.text)['status']
        print('Job status is {}'.format(status))

        if status == 'SUCCESS':
            ready = True
        elif status in ['FAILURE','REVOKED']:
            print('Error')
            break
        else:
            time.sleep(5)
    else:
        url_result = '{}/api/result/pdb_file/dssp/{}/'.format(dssp_url, job_id)

        response_result = requests.get(url_result)
        response_result.raise_for_status
        result = json.loads(response_result.text)['result']

        print(f"{structure_file.split('.')[0]}.dssp is created!")
        return result

def apply_dssp(structure_file):
    exit_code = os.system(f"mkdssp {structure_file} {structure_file.split('.del.pdb')[0]}.dssp")
    sys.exit(exit_code >> 8)

def main():
    uncorrected_dssp_files = []
    
    num = 1
    for structure in structure_list:
        structure_path = os.path.join(StructureSet_path, structure)
        preprocessing_path = os.path.join(structure_path, "preprocessing")
        
        structure_files = glob.glob(os.path.join(preprocessing_path, '*.del.pdb'))
        for structure_file in structure_files:
            structure_name = structure_file.split('/')[-1].split('.del.pdb')[0]
            #print(f"{num} --- {structure_name}")
            num += 1
            
            if args.locally:
                apply_dssp(structure_file)
                if not os.path.exists(os.path.join(preprocessing_path, f"{structure_name}.dssp")):
                    uncorrected_dssp_files.append(f"{structure_name}\tlocal problem")
            else:
                result = download_dssp(structure_file)    
                if result != None:
                    with open(os.path.join(preprocessing_path, f"{structure_name}.dssp"), 'w') as file:
                        file.write(result)
                else:
                    uncorrected_dssp_files.append(f"{structure_name}\tserver problem")

    if len(uncorrected_dssp_files) > 0:
        if os.path.exists(log_file):
            with open(log_file, 'a') as file:
                file.write('Not generated dssp files:\n' + '\n'.join(uncorrected_dssp_files) + '\n')
        else:
            with open(log_file, 'w') as file:
                file.write('Not generated dssp files\n' + '\n'.join(uncorrected_dssp_files) + '\n')

''' Launch script '''
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")
