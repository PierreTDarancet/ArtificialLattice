from types import SimpleNamespace
from ipywidgets import interact
import matplotlib
from matplotlib import pyplot
from mpl_toolkits import mplot3d
import numpy as np
import scipy.linalg as la 
import kwant
from kwant.wraparound import wraparound, plot_2d_bands

graphene = kwant.lattice.general([[np.sqrt(3)/2,1/2],[np.sqrt(3)/2,-1/2]],  #Lattice vectors 
                                  [[0,0],[1/np.sqrt(3),0]]) # Co-ordinates
a,b = graphene.sublattices 

Zigzag = kwant.lattice.general([[np.sqrt(3)/3,0],[0,1]], #Lattice vectors
                                     [[0,1/6],[np.sqrt(3)/2,2/6],[np.sqrt(3)/2,4/6],[0,5/6]])


def momentum_to_lattice(k):
    """Transform momentum to the basis of reciprocal lattice vectors.
    
    See https://en.wikipedia.org/wiki/Reciprocal_lattice#Generalization_of_a_dual_lattice
    """
    B = np.array(graphene.prim_vecs).T
    A = B.dot(np.linalg.inv(B.T.dot(B)))
    return np.linalg.solve(A, k)


def dispersion_2D(syst, args=None, lim=1.5*np.pi, num_points=200):
    """A simple plot of 2D band structure."""
    if args is None:
        args = []
    momenta = np.linspace(-lim, lim, num_points)
    energies = []
    for kx in momenta:
        for ky in momenta:
            lattice_k = momentum_to_lattice([kx, ky])
            h = syst.hamiltonian_submatrix(args=(list(args) + list(lattice_k)))
            energies.append(np.linalg.eigvalsh(h))
    
    energies = np.array(energies).reshape(num_points, num_points, -1)
    emin, emax = np.min(energies), np.max(energies)
    kx, ky = np.meshgrid(momenta, momenta)
    fig = pyplot.figure()
    axes = fig.add_subplot(1, 1, 1, projection='3d')
    for band in range(energies.shape[-1]):
        axes.plot_surface(kx, ky, energies[:, :, band], cstride=2, rstride=2,
                          cmap=matplotlib.cm.RdBu_r, vmin=emin, vmax=emax,
                          linewidth=0.1)

#Some helper functions to make the schematic look neat
def family_color(site): 
    if site.family == a: 
        return 0 
    else: 
        return 1

def same_color(site): 
    return 0
    
def hopping_lw(site1,site2): 
    return 0.1 if A in [site1.family,site2.family] else 0.05

def get_width(N):
    if N < 2:
        raise("N cannot be less than 2")
    else:
        return (N/4)*Zigzag.prim_vecs[1][1]

def get_length(L):
    if L < 2:
        raise("L cannot be less than 2")
    else:
        return (L/2)*Zigzag.prim_vecs[0][0]


def translate(pos,t,N,L):
    """ Translate the atoms within the unit cell keeping the boundries fixed
    Effectively changing the termination of the edges. Useful in scanning topological phases
    See S. Louie et. al PRL 2017; Mei-Yin Chou NanoLett 2018
    
    Parameters: 
    ========== 
        pos: position of the sites in the unit cell 
        N: number of atoms along the width 
        L: number of atoms along the length 
        
        N,L is required to map the atoms back into the unit cell when the are translated out of the cell
    TODO: 
    1. Make it more general for non-orthorhombic unit cells"""
    
    pos = np.array(pos)
    pos[:,0] += t 
    for p in pos: 
        if p[0] >= get_length(L): 
            p[0] -= get_length(L)
        elif p[0] <=0: 
            p[0] +=get_length(L) 
    return pos


def make_zigzag_ribbon(N=7, L = 5):
    """Returns a zigzag nanoribbon
    with length L and width of N carbon atoms

    Parameters:
    ==========

    N = Number of C atoms (width)
    L = Length of the nanoribbon along -x direction

    Returns:
    =======
    Instance of kwant.Builder() with the desired nanoribbon geometery
    """
    Zigzag = kwant.lattice.general([[np.sqrt(3)/3,0],[0,1]], #Lattice vectors
                                     [[0,1/6],[np.sqrt(3)/2,2/6],[np.sqrt(3)/2,4/6],[0,5/6]]) # Coordinates
    #Z_ribbon = kwant.Builder(kwant.TranslationalSymmetry([get_length(L),0]))

    Z_ribbon = kwant.Builder()
    Z_ribbon[Zigzag.shape((lambda pos: pos[1] >= 0 and pos[1] < get_width(N) and pos[0] >= 0 and pos[0] < get_length(L)), (0,0))] = 1
    Z_ribbon[Zigzag.neighbors()] = -1
    return Z_ribbon


def make_cove_edged_graphene(N=10,L=6):
    
    """Returns a zigzag cove edged graphene nanoribbon 
    with max width N carbon atoms and length of L carbon atoms
    
    Parameters:
    ==========
        N = Number of C atoms along the width including the edges 
        L = Number of C aroms along the length of the nanoribbon
    
    Returns:
    ========
    Instance of kwant.Builder() with cove edge geometry"""
    
    parent_ribbon = make_zigzag_ribbon(N=N,L=L)
    sites = list(parent_ribbon.sites())
    pos = []
    tag = []
    fam = []
    for site in sites:
        pos.append(site.pos)
        tag.append(site.tag)
        fam.append(site.family)
    pos = np.array(pos)
    # Get the 2nd smallest x values 
    x_values = np.unique(pos[:,0])
    y_values = np.unique(pos[:,1])
    
    max_x = x_values[x_values.argsort()[-1]]
    max2nd_x = x_values[x_values.argsort()[-2]]
    
    min_x = x_values[x_values.argsort()[0]]
    min2nd_x = x_values[x_values.argsort()[1]]
    
    min_y = y_values[y_values.argsort()[0]]
    min2nd_y = y_values[y_values.argsort()[1]]
    
    max_y = y_values[y_values.argsort()[-1]]
    max2nd_y = y_values[y_values.argsort()[-2]]
    
#    for p,t,f in zip(pos,tag,fam): 
#        if min2nd_x < p[0] < max_x: 
#            if abs(p[1]-min_y) < 1.e-3 or abs(p[1] - min2nd_y) < 1.0e-3:
#                del parent_ribbon[f(t[0],t[1])]
#        if min_x < p[0] < max2nd_x:
#            if abs(p[1]-max2nd_y) < 1.e-3 or abs(p[1] - max_y) < 1.0e-3:
#                del parent_ribbon[f(t[0],t[1])]
#    cove_ribbon = parent_ribbon
#    return cove_ribbon
    
    for p,t,f in zip(pos,tag,fam): 
        if  p[0] < x_values[x_values.argsort()[3]]:
            if abs(p[1]-min_y) < 1.e-3 or abs(p[1] - min2nd_y) < 1.0e-3:
                del parent_ribbon[f(t[0],t[1])]
        if min2nd_x < p[0] < max_x:
            if abs(p[1]-max2nd_y) < 1.e-3 or abs(p[1] - max_y) < 1.0e-3:
                del parent_ribbon[f(t[0],t[1])]
    cove_ribbon = parent_ribbon
    return cove_ribbon

def plot_wf(syst,i_start,i_end,ham):
    """Plot the wave function mapping on system with Hamiltonian
    "ham" in a PyWidget starting from band index i_start and
    ending at i_end"""
    eig_val,eig_vec = la.eigh(ham)
    def plot_band(i=0):
        print("Plotting wave function with index",i)
        print("Energy of the corresponding mode",eig_val[i], "x t")
        fig = kwant.plotter.map(syst,abs(eig_vec[:,i])**2,oversampling=50)
        fig.savefig('figures/%s.jpg'%i,dpi=400,quality=100,transparent=True)

    interact(plot_band,i=(i_start,i_end))

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
    system_1D[lattice_1D.neighbors()] = -1
    return system_1D 
