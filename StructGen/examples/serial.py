import kwant 
import numpy as np 
import math,os,datetime,sys 
import sys
import random
sys.path.append('/home/hyda/ArtificialLattice/StructGen/src')
from StructGen import StructGen,Check_redundant

N = int(sys.argv[1])
now = datetime.datetime.now() 
Struct_Dir='./RandomLineStructs-'+str(now).split()[0]+'/'
if not os.path.exists(Struct_Dir): 
    os.makedirs(Struct_Dir)

gen = StructGen(lx=6,ly=6) 
CR = Check_redundant()

trivial =[] 
nontrivial = []

def lines_sampling(): 
    ntrail_lines =0 
    while CR.is_redundant(gen.syst):
        ntrail_lines +=1
        gen.random_mirror_symmetric()
        if ntrail_lines > 1000: 
            swap_moves()

def swap_moves(): 
    ntrail_swap = 0 
    while CR.is_redundant(gen.syst): 
        ntrail_swap +=1 
        gen.swap_move()
        if ntrail_swap > 1000: 
            lines_sampling()


for n in range(N): 
    rand = random.uniform(0,1) 
    if rand < 0.5: 
        lines_sampling() 
    else: 
        swap_moves()
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

     
