# Generates NStruct syst (by rank=0) and calculates the Z2 in parallel

import sys
sys.path.append('/home/share/cnm50256/Quantum/ArtificialLattice/StructGen/src')
from StructGen import StructGen 
from mpi4py import MPI 
import os
from glob import glob
from math import floor

comm = MPI.COMM_WORLD 
rank = comm.Get_rank() 
nproc = comm.Get_size() 

NStruct = int(sys.argv[1])
nontrivial =[] 
trivial = []

gen = StructGen(lx=10,ly=10) 

#-------- Structure Generation----------#
if rank ==0:
    if not os.path.exists('./POSCARS'):
        os.makedirs('./POSCARS')
    for i in range(NStruct):
        status = gen.random_mirror_symmetric()
        if status:
            gen.syst2poscar('./POSCARS/POSCAR.{}'.format(i))
            syst = gen.get_syst()

#-------- Work load distributionn----------#
if rank ==0: 
    #Distribution
    structs_all = glob('./POSCARS/*')
    structs = [[] for _ in range(nproc)]
    j=0
    for struct in structs_all: 
        structs[j].append(struct) 
        j +=1
        if j > nproc-1: 
            j=0
    print(structs) 
    print(len(struct))
else: 
    structs=None
comm.Barrier() 
structs = comm.scatter(structs,root=0)
print(structs)

if rank ==0:
    if not os.path.exists('./Calc'):
        os.makedirs('./Calc')
out_file = open('./Calc/{}.out'.format(rank),'w')
for cand in structs:
    syst = gen.poscar2syst(cand) 
    index = int(cand.split('.')[-1])
    pol=gen.get_pol()
    Z2 = round(2*pol)%2
    if Z2 == 1: 
        nontrivial.append(cand)
        out_file.write("Struct index {} is non-trivial \n".format(index))
    elif Z2==0: 
        trivial.append(cand) 
        out_file.write("Struct index {} is trivial \n".format(index))
    else: 
        print("Z2 neither 1 or 0") 
        raise
    out_file.flush()
comm.Barrier() 

