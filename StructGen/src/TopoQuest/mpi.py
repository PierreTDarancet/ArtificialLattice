from StructGen import StructGen 
from mpi4py import MPI 
import sys
import os 

comm = MPI.COMM_WORLD 
rank = comm.Get_rank() 
nproc = comm.Get_size() 

#print(rank)

#NStruct = int(sys.argv[1])
##gen = StructGen(lx=25,ly=25) 
#gen = StructGen() 
#
#if rank ==0: 
#    candidates = []
#    for i in range(NStruct): 
#        gen.random_mirror_symmetric()
#        syst = gen.get_syst()
#        candidates.append(syst)
#        print(sys.getsizeof(syst))
#else: 
#    candidates=None 
#print(candidates)
##candidates = comm.scatter(candidates,root=0) 
#
##if rank ==0:
##    if not os.path.exists('./mpi'):
##        os.makedirs('./mpi')
#
#for cand in candidates: 
#    gen.set_syst(syst) 
#    pol=gen.get_pol()
#    print(rank,pol)

