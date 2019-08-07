from types import SimpleNamespace
from ipywidgets import interact
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import scipy.linalg as la 
from math import floor
import z2pack
from Z2_compute import zak_bands as zak_bands_z2
import kwant
from kwant.wraparound import wraparound, plot_2d_bands
from TopoQuest.StructGen import StructGen 

gen = StructGen('Armchair',nlx=1,nly=20)
Zigzag = kwant.lattice.general([[1,0],[0,np.sqrt(3)/3]], #Lattice vectors
                                     [[1/6,0],[2/6,np.sqrt(3)/2],[4/6,np.sqrt(3)/2],[5/6,0]]) # Coordinates

a,b,c,d = Zigzag.sublattices

def get_width(N=7): 
    if N < 2: 
        raise("N cannot be less than 2")
    else:
        return N/2*Zigzag.prim_vecs[1][1] + 0.01
    
def get_length(L=8): 
    return L/4*Zigzag.prim_vecs[0][0]
    
def make_1D_zigzag(N=7,L=8):
    #syst = kwant.Builder(kwant.TranslationalSymmetry(Zigzag.prim_vecs[0]))
    #syst = kwant.Builder(kwant.TranslationalSymmetry([get_length(L),0]))
    syst = kwant.Builder()
    syst[Zigzag.shape((lambda pos: pos[1] >0 and pos[1] <= get_width(N) and 0 <= pos[0] < get_length(L)),(0,0))] = 0
    syst[Zigzag.neighbors()] = -1
    return syst

def terminate_edges(syst): 
    sites = list(syst.sites())
    #print(sites)
    nsites = len(sites)
    pos = np.array([site.pos for site in sites])
    tags = [site.tag for site in sites]
    family = [site.family for site in sites]
    ymax = np.max(pos[:,1])
    ymin = np.min(pos[:,1])
    edge_index = []
    for i,p in enumerate(pos): 
        if abs(p[1] - ymax) < 1.e-2 or abs(p[1]-ymin) < 1.e-2: 
            edge_index.append(i)
    nedges = len(edge_index)
    #print(edge_index)
    edge_hoping_pairs = []
    for i in range(nedges): 
        site1 = sites[edge_index[i]]
        #print(site1.pos)
        neigh_sites = syst.neighbors(site1)
        for site2 in neigh_sites: 
            if abs(site2.pos[1] - ymax) < 1.0e-2  or abs(site2.pos[1]-ymin)< 1.0e-2:
                hop_pair = [site1.tag,site2.tag]
                pair_seen = hop_pair in edge_hoping_pairs
                if not pair_seen:
                    syst[site1,site2] = -1 - 0.5
                    edge_hoping_pairs.append([site1.tag,site2.tag])
                    #print(site1.pos)
                    #print(site2.pos)
    #print(edge_hoping_pairs)
    #print(edge_index)           
    return syst

def finite_to_1D(system,lat_vec,trans_sym_direction='x'): 
    """Adds a translational symmetry on a finite system
    Useful for making complex geometries eg: cove-edged and chevron ribbon
    
    Parameters
    ==========
        system: instance of the finite system
        lat_vec: lattice vector of the translational symmetry 
        trans_sym_direction: 'x' or 'y' , direction of the translational symmetry
        
    TODO: 
    1. Currently only works for orthorhombic unit cells 
    2. Get the onsite and hopping values directly from the passed system
        Currently hard set inside the code"""
    
    sites = list(system.sites())
    pos = [site.pos for site in sites] 
    if trans_sym_direction=='x':
        a = lat_vec
        b = max(np.array(pos)[:,1])
        trans_vec=[a,0]
    elif trans_sym_direction=='y':
        a = max(np.array(pos)[:,0])
        b = lat_vec
        trans_vec=[0,b]
    else: 
        raise #"Translation Symmetry direction should be 'x' or 'y'"
    lattice_1D = kwant.lattice.general([[a,0],[0,b]],pos)
    system_1D = kwant.Builder(kwant.TranslationalSymmetry(trans_vec))
    if trans_sym_direction=='x':
        system_1D[lattice_1D.shape((lambda pos: 0< pos[1] <= b),(0,0))]=0 
    if trans_sym_direction=='y':
        system_1D[lattice_1D.shape((lambda pos: 0< pos[0] <= a),(0,0))]=0 
    system_1D[lattice_1D.neighbors(2)] = -1
    return system_1D 

def get_pol_custom(syst):
    pos = [ site.pos for  site in syst.sites()]
    nsites = len(pos)
    syst = syst.finalized()
    kwant.plotter.bands(syst)
    lattice = kwant.lattice.general([[get_length(4),0],[0,get_width(16)]],pos)
    act_pos = np.array([syst.pos(i) for i in range(nsites)])
    a1,a2 = [lattice.prim_vecs[0][0],get_width(N=7)]
    red_pos = np.zeros(np.shape(act_pos))
    red_pos[:,0] = act_pos[:,0]/a1
    red_pos[:,1] = act_pos[:,1]/a2


    ham_k=zak_bands_z2(syst,momenta=1001,dim=2)
    z2_system = z2pack.hm.System(ham_k,dim=2,#pos=red_pos,
                                     convention=1)
    result = z2pack.line.run(system=z2_system, 
                            line=lambda t1: [t1,0])#,#,n/(N+1)])#,
                            #pos_tol=1e-3,iterator=range(200,500,2));

    print("Polarization:",result.pol)
    
def make_junction(syst1,syst2,lx1,lx2,ly,xoff=0,yoff=0):
    lx = lx1 + lx2
    nx = 1
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
    lat = kwant.lattice.general([[nx*lx,0],[0,2*ly]],pos_all,norbs=1)
    syst = kwant.Builder(kwant.TranslationalSymmetry([lx,0]))
    #syst = kwant.Builder()
    syst[lat.shape((lambda pos: 0< pos[0]< nx*lx and min_y - yoff<=pos[1]<ly+yoff),(0,0))] = 0 
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
        fig = kwant.plotter.map(syst,abs(eig_vec[:,i])**2,oversampling=3,fig_size=(16,10))
        
        #fig.savefig('figures/%s.jpg'%i,dpi=400,quality=100,transparent=True)
    
    interact(plot_band,i=(i_start,i_end))
    
def check_junction(syst1,syst2,lx1,lx2,ly=14,xoff=0,yoff=0):
    syst = make_junction(syst1,syst2,lx1,lx2,ly,xoff,yoff)
    #syst.eradicate_dangling()
    kwant.plot(syst,site_color='black',fig_size=(16,10));
    
    nbands = len(syst.sites())
    n1 = int(nbands/2 -5)
    n2 = int(nbands/2 +5)
    print('Number of bands in junction is {}'.format(nbands))
    
    syst = syst.finalized()
    ham = syst.hamiltonian_submatrix() 
    fig = plot_wf(syst,n1,n2,ham)

wide = make_1D_zigzag(N=9,L=4)
wide = finite_to_1D(wide,1)
gen.syst = wide
gen.syst2poscar('POSCAR.wide')

narrow = make_1D_zigzag(N=7,L=4)
narrow = finite_to_1D(narrow,1)
gen.syst = narrow
gen.syst2poscar('POSCAR.narrow')


momenta = np.linspace(-np.pi,np.pi,101) 
bands1 = kwant.physics.Bands(wide.finalized())
energies1 = [bands1(k) for k in momenta] 
bands2 = kwant.physics.Bands(narrow.finalized())
energies2 = [bands2(k) for k in momenta] 

fig = plt.figure()
ax = fig.add_subplot(111) 

ax.plot(momenta,energies1,color='black',label='Wide')
ax.plot(momenta,energies2,color='red',label='Narrow')
plt.xticks(ticks=[-np.pi,0.0,np.pi],labels=[r'-$\pi$/a','0',r'$\pi$/a']) 
plt.yticks(ticks=[-1.0,-0.5,0.0,0.5,1.0])
    
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(20)
    #tick.label.set_fontsize(20)

for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(20)

plt.ylabel('Energy (x t) [eV]',fontweight='bold',fontsize=25)
plt.ylim([-1.0,1.0])
plt.xlim([-np.pi,np.pi])
#plt.legend()
plt.tight_layout()
#plt.show()
plt.savefig('bands.tiff')
fig.clear()                               

#gen2 = StructGen('Armchair',nlx=2,nly=20) 
#comp_syst = make_junction(wide,narrow,1,1,12,yoff=4*np.sqrt(3)/3)
#gen2.syst = comp_syst
#gen2.syst2poscar('POSCAR.junction')
#bands = kwant.physics.Bands(comp_syst.finalized()) 
#energies = [bands(k) for k in momenta]

def make_junction_comp(systs,lx1,yoffs):
    lx = len(systs)*lx1
    pos_all = []
    for i,syst in enumerate(systs):
        pos =  [[site.pos[0] + i*lx1,site.pos[1]+yoffs[i]] for site in syst.sites()]
        pos_all += pos
    min_y = np.min(np.array(pos_all)[:,1])
    max_y = np.max(np.array(pos_all)[:,1])
    lat = kwant.lattice.general([[lx,0],[0,max_y+5.0]],pos_all,norbs=1)
    syst = kwant.Builder(kwant.TranslationalSymmetry([lx,0]))
    syst[lat.shape((lambda pos: 0< pos[0]< lx and 0<=pos[1]<max_y+5.0),(0,0))] = 0 
    syst[lat.neighbors()] = -1 
    return syst 

systs1 = ([narrow] + [wide])*4 + [narrow] 
#systs = [narrow,narrow,wide,narrow,wide,narrow]
#yoffs = np.array([0,2,3,6,6,5,2,0])*np.sqrt(3)/3 
yoffs = np.array([0,1,2,3,4,3,2,1,0])*np.sqrt(3)/3 
comp_syst1 = make_junction_comp(systs1,1,yoffs)
#kwant.plot(comp_syst1,site_color='black')
momenta2 = np.linspace(-np.pi,np.pi,101) 
bands3 = kwant.physics.Bands(comp_syst1.finalized()) 
energies3 = [bands3(k) for k in momenta]

gen2 = StructGen('Armchair',nlx=2,nly=20) 
gen2.syst = comp_syst1
gen2.syst2poscar('POSCAR.wiggle')


systs2 = [wide]*9
#systs = [narrow,narrow,wide,narrow,wide,narrow]
#yoffs = np.array([0,2,3,6,6,5,2,0])*np.sqrt(3)/3 
yoffs = np.array([0]*12)*np.sqrt(3)/3 
comp_syst2 = make_junction_comp(systs2,1,yoffs)
#kwant.plot(comp_syst2,site_color='black')
bands4 = kwant.physics.Bands(comp_syst2.finalized()) 
energies4 = [bands4(k) for k in momenta2]

systs3 = [narrow]*9
#systs = [narrow,narrow,wide,narrow,wide,narrow]
#yoffs = np.array([0,2,3,6,6,5,2,0])*np.sqrt(3)/3 
yoffs = np.array([0]*12)*np.sqrt(3)/3 
comp_syst3 = make_junction_comp(systs3,1,yoffs)
#kwant.plot(comp_syst3,site_color='black')
bands5 = kwant.physics.Bands(comp_syst3.finalized()) 
energies5 = [bands5(k) for k in momenta2]

fig = plt.figure()
ax = fig.add_subplot(111) 

ax.plot(momenta2,energies3,color='black',label='Wide')
ax.plot(momenta2,energies4,color='red',label='Wide')
ax.plot(momenta2,energies5,color='blue',label='Wide')
plt.xticks(ticks=[-np.pi,0.0,np.pi],labels=[r'-$\pi$/a','0',r'$\pi$/a']) 
plt.yticks(ticks=[-1.0,-0.5,0.0,0.5,1.0])
    
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(20)
    #tick.label.set_fontsize(20)

for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(20)

plt.ylabel('Energy (x t) [eV]',fontweight='bold',fontsize=25)
plt.ylim([-1.0,1.0])
plt.xlim([-np.pi,np.pi])
#plt.legend()
plt.tight_layout()
plt.savefig('bands-comp-syst.tiff')
fig.clear()                               

#gen2.syst = comp_syst
#gen2.syst2poscar('POSCAR.junction')
#bands = kwant.physics.Bands(comp_syst.finalized())

