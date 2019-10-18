import sys 
sys.path.append('/Users/ssrinivasan/ArtificialLattice/StructGen/src/')
sys.path.append('/Users/ssrinivasan/transfer/4x4/sandbox/structs/4x4/')
from glob import glob
import numpy as np 
from TopoQuest.utilities import terminate_edges
from TopoQuest.StructGen import StructGen 
import random
import kwant
import scipy.linalg as la 
import matplotlib.pyplot as plt


def make_junction(syst1,syst2,lx1,lx2,ly,xoff=0,yoff=0):
    lx = lx1 + lx2
    nx = 3
    nx1 = nx 
    nx2 = nx
    pos_all0 = []
    for i in range(nx1):
        pos1 =  [[site.pos[0] + i*lx1,site.pos[1]] for site in syst1.sites()]
        pos_all0 += pos1
    for i in range(nx2): 
        pos2 =  [[site.pos[0] + nx1*lx1 + i*lx2+xoff,site.pos[1]+yoff] for site in syst2.sites()]
        pos_all0 += pos2
    pos_all = pos_all0

    min_y = np.min(np.array(pos_all)[:,1])
    lat = kwant.lattice.general([[nx*lx,0],[0,3*ly]],pos_all,norbs=1)
    syst = kwant.Builder()
    syst[lat.shape((lambda pos: 0< pos[0]< nx*lx+xoff and min_y <=pos[1]<ly),(0,0))] = 0 
    syst[lat.neighbors()] = -1 
    print(len(pos1),len(pos2),len(pos_all))
    return syst 

def plot_wf(syst,i_start,i_end,ham):
    """Plot the wave function mapping on system with Hamiltonian 
    "ham" in a PyWidget starting from band index i_start and 
    ending at i_end"""
    eig_val,eig_vec = la.eigh(ham)
    def plot_band(i=0): 
        print("Plotting wave function with index",i)
        print("Energy of the corresponding mode",eig_val[i], "x t")
        fig = plt.figure(figsize=(16,3))
        kwant.plotter.map(syst,abs(eig_vec[:,i])**2,oversampling=3,fig_size=(16,10))
        plt.axis('off')
        fig.savefig('figures/%s.jpg'%i,dpi=600,quality=300,transparent=True)
    interact(plot_band,i=(i_start,i_end))
    
def check_junction(syst1,syst2,lx1,lx2,ly=14,xoff=0,yoff=0):
    syst = make_junction(syst1,syst2,lx1,lx2,ly,xoff,yoff)
    #syst.eradicate_dangling()
    #kwant.plot(syst,site_color='black');
    
    nbands = len(syst.sites())
    n1 = int(nbands/2 -5)
    n2 = int(nbands/2 +5)
    print('Number of bands in junction is {}'.format(nbands))
    #syst_fin = syst.finalized()
    #ham = syst_fin.hamiltonian_submatrix() 
    #fig = plot_wf(syst_fin,n1,n2,ham)
    return syst


base = '/Users/ssrinivasan/transfer/4x4/sandbox/structs/4x4/'
nontrivial = np.loadtxt('/Users/ssrinivasan/Quantum/4x4/sandbox/orig_z2.dat',usecols=1)
full_id_list = np.array([i for i in range(100)])
trivial = [item for item in full_id_list if item not in nontrivial]
rnt = base+'POSCAR.37'
rt = base+'POSCAR.28'

gen_nt = StructGen('Armchair',nlx=4,nly=4) 
gen_nt.poscar2syst(rnt)
nt_syst = gen_nt.get_syst()

gen_tt = StructGen('Armchair',nlx=4,nly=4)
gen_tt.poscar2syst(rt)
tt_syst = gen_tt.get_syst()

tt_syst_pos = np.array([site.pos for site in tt_syst.sites()])
yoff = np.min(tt_syst_pos[:,1])
syst = check_junction(nt_syst,tt_syst,lx1=4,lx2=4,ly=4,yoff=-1*yoff)

gen_tt.syst = syst 
gen_tt.syst2poscar(filename='POSCAR.junc')
