#import kwant 
#import numpy as np 
import os,datetime,sys 
import networkx as nx 
#import math
#import sys
#import random
from TopoQuest.StructGen import StructGen,Check_redundant

N = int(sys.argv[1])
now = datetime.datetime.now() 
Struct_Dir='./swap-'+str(now).split()[0]+'/'
if not os.path.exists(Struct_Dir): 
    os.makedirs(Struct_Dir)

gen = StructGen('Armchair',nlx=4,nly=8,dumpfile=Struct_Dir+'dump.atom') 
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
            print('lines sampling')
            lines_sampling()


for n in range(N):
    try:
        swap_moves()
    except: 
        #pos_gen = gen._get_site_pos()
        #pos_double = gen._get_site_pos(syst=gen.full_double_syst)
        #nx.draw_networkx(gen.graph,pos=pos_gen)    
        #nx.draw_networkx(gen.full_double_graph,pos=pos_double)
        raise
    gen.syst2dump(frame=n)


     
