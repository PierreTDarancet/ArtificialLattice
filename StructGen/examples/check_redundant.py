import sys,os 
from shutil import copyfile 
from glob import glob 
sys.path.append('/Users/ssrinivasan/ArtificialLattice/StructGen/src')
from StructGen import StructGen, Check_redundant 

poscar_dirs = sys.argv[1:]
POSCARS = []
for direc in poscar_dirs: 
    for item in glob(direc+'/*'): 
        POSCARS.append(item)

outdir = './Unique_structs/' 
if not os.path.exists(outdir):
    os.makedirs(outdir)

gen = StructGen(lx=6,ly=6) #Needs to match with the POSCAR files being read 
CR = Check_redundant()
purge_list = []
exceptions = []
for POSCAR in POSCARS:
    print(POSCAR)
    try:
        syst = gen.poscar2syst(POSCAR)
        is_redundant = CR.is_redundant(gen) 
        if is_redundant:
            purge_list.append(POSCAR) 
    except: 
        exceptions.append(POSCAR)
        #raise
n_unqiue_structs = len(POSCARS) - len(purge_list)
cnt = 0 
for POSCAR in POSCARS: 
    if POSCAR not in purge_list:
        copyfile(POSCAR,outdir+'POSCAR.{}'.format(cnt))
        cnt +=1 

with open('exceptions.dat','w') as f: 
    for item in exceptions: 
        f.write(item+'\n') 

        


