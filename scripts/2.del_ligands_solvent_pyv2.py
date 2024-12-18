from chimera import runCommand as run
import os
import glob

''' Set paths '''
main_path = os.getcwd()
StructureSet_path = os.path.join(main_path, "results")
structure_list = [i for i in os.listdir(StructureSet_path) if '.' not in i]

c = 0
for structure in sorted(structure_list):  
    structure_path = os.path.join(StructureSet_path, structure)
    preprocessing_path = os.path.join(structure_path, "preprocessing")
    structure_files = glob.glob(os.path.join(preprocessing_path, "*.pdb"))
    
    for structure_file in structure_files:
        c += 1
        #print "{0} --- {1}".format(c, structure_file.split('.')[0].split('/')[-1])
        
        ''' Chimera Run '''
        run("open {}".format(structure_file))
        run("del ligand")
        run("del solvent")
        run("write #0 {}".format(structure_file))
        run("del")

        ''' Processing of output PDB file from chimera: two space must be added! '''
        with open(structure_file, 'r') as file:
            data = file.read()

        new_ATOM = []
        for stroka in data.strip('\n').split('\n'):
            if "ATOM" in stroka or "HETATM" in stroka:
                new_stroka = stroka + ' '*2
                new_ATOM.append(new_stroka)

        with open(structure_file.split('.pdb')[0] + '.del.pdb', 'w') as file:
            file.write('\n'.join(new_ATOM))