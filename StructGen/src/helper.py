#from ipywidgets import interact
import matplotlib
from matplotlib import pyplot
import numpy as np
import scipy.linalg as la 
import kwant
import functools

graphene = kwant.lattice.general([[np.sqrt(3)/2,1/2],[np.sqrt(3)/2,-1/2]],  #Lattice vectors 
                                  [[0,0],[1/np.sqrt(3),0]]) # Co-ordinates
a,b = graphene.sublattices 

Zigzag = kwant.lattice.general([[np.sqrt(3)/3,0],[0,1]], #Lattice vectors
                                     [[0,1/6],[np.sqrt(3)/2,2/6],[np.sqrt(3)/2,4/6],[0,5/6]])

Armchair = kwant.lattice.general([[1,0],[0,np.sqrt(3)/3.0]], #Lattice vectors
                                     [[1/6,0],[2/6,np.sqrt(3)/2],[4/6,np.sqrt(3)/2],[5/6,0]])

Armchair_trans = kwant.lattice.general([[1,0],[0,np.sqrt(3)/3.0]], #Lattice vectors
                                    [[1/6,0],[2/6,-1*np.sqrt(3)/2],[4/6,-1*np.sqrt(3)/2],[5/6,0]])


lat_dict = {'armchair':Armchair,'zigzag':Zigzag,
            'armchair_trans':Armchair_trans}
trans_vec_dict = {'armchair':Armchair.prim_vecs[0],'zigzag':Zigzag.prim_vecs[1],
                  'armchair_trans':Armchair_trans.prim_vecs[0]}

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

   
def hopping_lw(syst):
    def hopping_lw_by_overlap(site1,site2,syst):    
        abshop = abs(syst[site1,site2])
        if abshop <=1: 
            return 0.1 
        else:
            return 0.15
    return functools.partial(hopping_lw_by_overlap,syst=syst)
    

def get_width(N,lat):
    num_c_atoms_per_cell = {Armchair:2,Zigzag:4,
                            Armchair_trans:2} 

    return float((N/num_c_atoms_per_cell[lat]))*lat.prim_vecs[1][1]

def get_length(L,lat):
    num_c_atoms_per_cell = {Armchair:4,Zigzag:2,
                            Armchair_trans:4}
    if L < 2:
        raise("L cannot be less than 2")
    else:
        return (L/num_c_atoms_per_cell[lat])*lat.prim_vecs[0][0]


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


def make_1D_cell(N=7,edge='armchair',offset=0):
    lat = lat_dict[edge] 
    trans_vec = trans_vec_dict[edge]
    syst = kwant.Builder(kwant.TranslationalSymmetry(trans_vec))
    #syst= kwant.Builder()
    syst[lat.shape((lambda pos: pos[1] >= offset and 
                    pos[1] < get_width(N,lat)+offset),(0,0))] = 0
    syst[lat.neighbors()] = -1
    return syst


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
    site_pos = [site.pos for site in sites] 
    if trans_sym_direction=='x':
        a = lat_vec
        b = max(np.array(site_pos)[:,1])
        trans_vec=[a,0]
    elif trans_sym_direction=='y':
        a = max(np.array(site_pos)[:,0])
        b = lat_vec
        trans_vec=[0,b]
    else: 
        raise #"Translation Symmetry direction should be 'x' or 'y'"
    lattice_1D = kwant.lattice.general([[a,0],[0,b]],site_pos)
    system_1D = kwant.Builder(kwant.TranslationalSymmetry(trans_vec))
    
    def fill(pos): 
        if pos in site_pos: 
            return True 
        else: 
            return False
        
    system_1D[lattice_1D.shape(fill,(0,0))]=0 
    system_1D[lattice_1D.neighbors()] = -1
    return system_1D 

def terminate_edges(syst,delta_t, num_bulk_neigh=3): 
    """ The overlap or hopping between the edge sties is increased by "delta_t" 
    effectively performing termination of the sites at the edges to match the 
    site coodrination of bulk sites
    
    Parameters:
    ==========
        syst: kwant.Builder instance
        delta_t: change in overlap after termination 
        num_bulk_neigh: number of neighbors for bulk sites
    """
    sites = list(syst.sites())
    site_degree = [syst.degree(site) for site in sites]
    edge_sites=[]
    site_seen = []
    for i,site in enumerate(sites): 
        if site_degree[i] < num_bulk_neigh:
            edge_sites.append(site)
    for edge_site in edge_sites:  
        neighbors = syst.neighbors(edge_site)
        for neigh in neighbors: 
            if syst.degree(neigh) < num_bulk_neigh:
                if neigh.tag not in site_seen:
                    # Hopping for the same pair is increased twice 
                    # so we only increase half each time 
                    syst[edge_site,neigh] += delta_t/2.0  
        site_seen.append(edge_site.tag)
    return syst



