#import kwant 
#import numpy as np 
import os,datetime,sys 
import networkx as nx 
#import math
#import sys
#import random
from TopoQuest.StructGen import StructGen,Check_redundant

N = 500
now = datetime.datetime.now() 
Struct_Dir='./swap-'+str(now).split()[0]+'/'
if not os.path.exists(Struct_Dir): 
    os.makedirs(Struct_Dir)

gen = StructGen('Armchair',nlx=4,nly=4,dumpfile=Struct_Dir+'dump.atom') 
CR = Check_redundant()
gen.random_inversion_symmetric()

trivial =[] 
nontrivial = []

def lines_sampling(): 
    ntrail_lines =0 
    while CR.is_redundant(gen.syst):
        ntrail_lines +=1
        gen.random_inversion_symmetric()
        if ntrail_lines > 1000: 
            swap_moves(sym='inversion')

def swap_moves(): 
    ntrail_swap = 0
    ntrail_lines = 0
    gen.swap_move(sym='inversion')
    while CR.is_redundant(gen.syst): 
        ntrail_swap +=1 
        if ntrail_swap > 1000: 
            ntrail_lines +=1
            if ntrail_lines > 100: 
                gen._fill_all_sites()
            lines_sampling()
        gen.swap_move(sym='inversion')



for n in range(N):
    try:
        swap_moves()
    except: 
        #pos_gen = gen._get_site_pos()
        #pos_double = gen._get_site_pos(syst=gen.full_double_syst)
        #nx.draw_networkx(gen.graph,pos=pos_gen)    
        #nx.draw_networkx(gen.full_double_graph,pos=pos_double)
        #gen.plot_syst()
        #gen.draw_lattice_graph()
        raise
    gen.syst2dump(frame=n)


     
