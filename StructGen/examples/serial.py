from StructGen import StructGen
import kwant 
import numpy as np 
import math,os,datetime,sys 

N = int(sys.argv[1])
width = 5
now = datetime.datetime.now() 
Struct_Dir='./RandomLineStructs-'+str(now).split()[0]+'/'
if not os.path.exists(Struct_Dir): 
    os.makedirs(Struct_Dir)

gen = StructGen() 

trivial =[] 
nontrivial = []
for n in range(N): 
    gen.random_mirror_symmetric(Ncentral=5) 
    gen.syst2poscar(Struct_Dir+'POSCAR.'+str(n))
    pol = gen.get_pol() 
    Z2 = round(2*pol)%2
    if Z2 == 1: 
       nontrivial.append(n)
    elif Z2 == 0: 
       trivial.append(n)
    else: 
       print("Z2 neither 1 or 0")
       raise
with open("nontrivial.dat","w") as f: 
    f.write("\t".join([str(item) for item in nontrivial]))

     
