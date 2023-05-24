''' Import libraries '''
import os
import glob
import time
import requests
import json
import argparse

''' Terminal '''
parser = argparse.ArgumentParser()
parser.add_argument("-l", "--locally", help="Run DSSP locally (specify '-l') or via server", action="store_true")
args = parser.parse_args()

''' Set paths '''
main_path = os.getcwd()
StructureSet_path = os.path.join(main_path, "structures")
structure_list = [i for i in os.listdir(StructureSet_path) if '.' not in i]

def download_dssp(structure_file):
    dssp_url = 'https://www3.cmbi.umcn.nl/xssp/'
    url_create = '{}api/create/pdb_file/dssp/'.format(dssp_url)

    response_create = requests.post(url_create, files={'file_':open(structure_file, 'rb')})
    response_create.raise_for_status

    job_id = json.loads(response_create.text)['id']
    print("Job submitted successfully. Id is: '{}'".format(job_id))

    ready = False
    while not ready:

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
    os.system(f"mkdssp -i {structure_file} -o {structure_file.split('.pdb')[0]}.dssp")

def main():
    uncorrected_dssp_files = []
    
    num = 1
    for structure in structure_list:
        structure_path = os.path.join(StructureSet_path, structure)
        
        structure_files = glob.glob(os.path.join(structure_path, '*_*.pdb')) + glob.glob(os.path.join(structure_path, "AF-*-F1.pdb"))
        for structure_file in structure_files:
            structure_name = structure_file.split('/')[-1].split('.')[0]
            print(f"{num} --- {structure_name}")
            num += 1
            
            if args.locally:
                apply_dssp(structure_file)
                if not os.path.exists(f"{structure_file.split('.pdb')[0]}.dssp"):
                    uncorrected_dssp_files.append(f"{structure_name}.pdb - error! (locally)")
            else:
                result = download_dssp(structure_file)
                if result != None:
                    with open(f"{structure_file.split('.')[0]}.dssp", 'w') as file:
                        file.write(result)
                else:
                    uncorrected_dssp_files.append(f"{structure_name}.pdb --- error! (via server)")

    if len(uncorrected_dssp_files) > 0:
        with open(os.path.join(StructureSet_path, f"UncorrectedDSSPFiles.txt"), 'w') as file:
            file.write('\n'.join(uncorrected_dssp_files))

''' Launch script '''
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")
