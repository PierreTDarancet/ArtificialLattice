#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 16:27:19 2019

@author: ssrinivasan
"""
#from helper import *
from .utilities import Armchair,get_width
from pymatgen.io.lammps.data import LammpsData 
from pymatgen import Lattice, Structure
from math import floor,ceil,copysign
from operator import itemgetter
import kwant,z2pack
import cmath,random,operator 
import numpy as np 
import networkx as nx 
import networkx.algorithms.isomorphism as iso
import itertools
import copy

def get_Hk(sys, args=(), momenta=65, file=None, *, params=None,dim=3):
    """Returns hamiltonian as a function of k. Modified from kwant's 
    plot.spectrum

    Parameters
    ----------
    sys : kwant.system.InfiniteSystem
        A system bands of which are to be plotted.
    args : tuple, defaults to empty
        Positional arguments to pass to the ``hamiltonian`` method.
        Mutally exclusive with 'params'.
    momenta : int or 1D array-like
        Either a number of sampling points on the interval [-pi, pi], or an
        array of points at which the band structure has to be evaluated.
    file : string or file object or `None`
        The output file.  If `None`, output will be shown instead.
    params : dict, optional
        Dictionary of parameter names and their values. Mutually exclusive
        w
        ith 'args'.
    dim: int , optional
        Dimension of the system
        Default = 3 

    Returns
    -------
        h_k: A callable function defition which returns the hamiltonian 
            as a function of 'k'
    """


    syst = sys  # for naming consistency inside function bodies
    momenta = np.array(momenta)
    if momenta.ndim != 1:
        momenta = np.linspace(-np.pi, np.pi, momenta)

    # expand out the contents of 'physics.Bands' to get the H(k),
    # because 'spectrum' already does the diagonalisation.
    ham = syst.cell_hamiltonian(args, params=params)
    if not np.allclose(ham, ham.conjugate().transpose()):
        raise ValueError('The cell Hamiltonian is not Hermitian.')
    _hop = syst.inter_cell_hopping(args, params=params)
    hop = np.empty(ham.shape, dtype=complex)
    hop[:, :_hop.shape[1]] = _hop
    hop[:, _hop.shape[1]:] = 0

    def h_k(k):
        # H_k = H_0 + V e^-ik + V^\dagger e^ik
        mat = hop * cmath.exp(-1j * np.array(2*np.pi*k[0]))
        mat +=  mat.conjugate().transpose() + ham
        return mat

    return h_k


def calc_pol(syst,red_pos,wcc=False ):   
    """ Returns the polarization of a unit cell built using kwant
    
    Parameters
    ----------
    syst: finalized kwant builder 
    red_pos: reduced position of the sites in the unit cell
    
    Returns 
    -------
    result.pol: Polarization from the result of the line calculation using Z2pack
    """
    
    Hk = get_Hk(syst,dim=2)
    z2_system = z2pack.hm.System(Hk,dim=2,#pos=red_pos,
                                 convention=2)
    result = z2pack.line.run(system=z2_system, 
                              line=lambda t1: [t1,0],
                              pos_tol=1e-2,iterator=range(7,501,2));
    if wcc: 
        return result.wcc, result.pol
    else: 
        return result.pol
    
class Check_redundant(): 
    """ Class containing methods to check for redundancy between the randomly
    generated structures by StructGen()
    
    Attributes 
    ---------- 
    seen_lat: set of set of sites of previously seen unit cells  
    """
    def __init__(self):
        self.seen_lat = set()
    
    def is_redundant(self,syst):
        """Returns True if the structure has been seen before by the generator
        
        Parameters
        ----------
        syst: Predicted builder by the generator 
        
        Returns 
        -------
        Boolean 
        """
        
        if syst is None: 
            return True 
        else:
            sites = frozenset(site for site in list(syst.sites()))
            if sites in self.seen_lat: 
                print("Redundant structure identified")
                return True 
            else: 
                self.seen_lat.add(sites)
                return False

class StructGen(): 
    """
    Class for random 1D structure generator 
    
    Attributes 
    ----------
    lat: kwant.lattice 
        lattice used for structure generation 
    lx: float 
        lenght of the unit cell 
    ly: float 
        width of the unit cell 
    onsite: float 
            onsite energy 
    hop: float 
         hopping 
    syst: kwant.Builder 
        geometry generated by the Class
    """
    
    def __init__(self,lat=Armchair,lx=5*Armchair.prim_vecs[0][0],
                 ly= 5*Armchair.prim_vecs[1][1],onsite=0,hop=-1,syst=None): 
        self.lat = lat
        self.lx  = lx
        self.ly  = ly
        self.onsite = onsite
        self.hop = hop
        self.syst = syst
        self.full_syst_pos = None 
        self.index_dict = None 
        self.full_graph = None
        self.full_double_graph = None 
        self.full_syst = None
        self.full_double_syst = None 
        self.mirror_sym_pairs = None 
        self._set_full_syst_attributes()
        self.graph = None 
        
    def _get_mirror_symmetric_pairs(self): 
        mirror_plane = self.lx/2.0
        mirror_sym_pairs = {}
        for site1 in list(self.full_double_syst.sites()): 
            x1,y1 = site1.pos
            for site2 in list(self.full_double_syst.sites()): 
                if site1 is not site2: 
                    x2,y2 = site2.pos
                    x_check = abs(abs(x1-mirror_plane)-abs(x2-mirror_plane)) < 1.e-2 
                    y_check = abs(y1-y2) < 1.e-2
                    if x_check and y_check:
                        mirror_sym_pairs[site1] = site2
        return mirror_sym_pairs
                                                
    def _get_site_pos(self,syst=None):
        if syst is None:
            pos = np.array([site.pos for site in list(self.syst.sites())]).tolist()
            return pos
        else: 
            pos = np.array([site.pos for site in list(syst.sites())]).tolist()
            return pos
        
    def _set_full_syst_attributes(self): 
        full_syst = self._make_full_syst(uly=0)
        pos = self._get_site_pos(syst=full_syst)
        pos.sort(key=operator.itemgetter(0,1))
        pos = np.array(pos)
        index_dict = {}
        for i,site in enumerate(pos): 
            index_dict[tuple(site)] =i
        self.syst = full_syst
        self.full_syst_pos = pos 
        self.index_dict = index_dict
        self.full_graph = self._construct_full_graph()
        self.full_syst = full_syst
        del full_syst 
        self.syst = None
        double_full_syst = self._make_full_syst(uly=-1*self.ly)
        self.syst = double_full_syst
        self.full_double_syst = double_full_syst
        self.full_double_graph = self._construct_full_graph() 
        del double_full_syst
        self.syst = None
        self.mirror_sym_pairs = self._get_mirror_symmetric_pairs()
        
    def _make_full_syst(self,uly): 
        full_syst = kwant.Builder()
        full_syst[self.lat.shape(
                (lambda pos: 0<=pos[0]<=self.lx and uly<=pos[1]<=self.ly),
                (0,0))] = self.onsite 
        full_syst[self.lat.neighbors()] = self.hop
        return full_syst
    
    def set_cell_attributes(self,lx=5,ly=5,lat=Armchair): 
        self.lx = lx 
        self.ly = ly 
        self.lat = lat
        self._set_full_syst_attributes() 
    
    def fill_all_sites(self): 
        full_syst=self._make_full_syst() 
        self.syst = full_syst
        
    def from_lmpdata(self,data):


        lmpdata = LammpsData.from_file(data)
        struct = lmpdata.structure
        topo = lmpdata.topology
        atoms = lmpdata.atoms
        bonds = topo['Bonds']
        nbonds = len(bonds)
        [a,b,c] = struct.lattice._matrix
        vec1 = [a[0],a[1]]
        vec2 = [b[0]*2,b[1]*2]
        pos = []

        
        for i in range(1,len(atoms)+1): 
            x = atoms.loc[i,'x']
            y = atoms.loc[i,'y']
            pos.append([x,y])
        pos = np.array(pos)
        max_y = np.max(pos[:,1])
        lat = kwant.lattice.general([vec1,vec2],pos)
        syst = kwant.Builder(kwant.TranslationalSymmetry(vec1))
        syst[lat.shape((lambda pos: 0<=pos[1]<= max_y),(0,0))]=self.onsite
        self.lat = lat
        self.lx = a[0]
        self.ly = max_y+0.5
        
        for i in range(1,nbonds+1): 
            at1 = atoms.loc[bonds.loc[i,'atom1']]
            at2 = atoms.loc[bonds.loc[i,'atom2']]
            x1,y1 = at1['x'],at1['y']
            x2,y2 = at2['x'],at2['y']
            site1 = syst.closest([x1,y1])
            site2 = syst.closest([x2,y2])
            syst[site1,site2] = self.hop
        syst[lat.neighbors()] = self.hop
        syst[lat.neighbors(2)] = self.hop
        syst[lat.neighbors(3)] = self.hop
        self.syst = syst 
        return syst 
        #min_x = np.min(pos[:,0])
        #max_x = np.min(pos[:,0])
        
    def _construct_full_graph(self,draw=False): 
              
        G = nx.Graph()
        for hopping,value in self.syst.hopping_value_pairs():
            u,v = hopping
            G.add_node(u,x=u.pos[0],y=u.pos[1]) 
            G.add_node(v,x=v.pos[0],y=u.pos[1]) 
            G.add_edge(u,v)
        return G 
        
        

        

    
    def _is_continous(self): 
        G = self.construct_graph() 
        graph_xlist = [node for node in G.nodes('x')]
        graph_xlist = sorted(graph_xlist,key=itemgetter(1))
        head,tail = graph_xlist[0][0],graph_xlist[-1][0]
        return nx.has_path(G,head,tail)
    
    def swap_move(self):
        
        
        edge_sites = []
        for site in self.syst.sites(): 
            if self.syst.degree(site)<3:
                edge_sites.append(site)
        possible_head_tails = []
        for (s,t) in itertools.combinations(edge_sites,2):
            if nx.shortest_path_length(self.full_double_graph,s,t) <7:
                possible_head_tails.append([s,t])
                

        def _add_ring(paths):
            pos = np.array(self._get_site_pos())
            ymin,ymax = np.min(pos[:,1]),np.max(pos[:,1])
            path_not_exist = []
            for path in paths: 
                if len(path) <7: 
                    site_exist_in_path = True 
                    for site in path: 
                        if site not in self.syst.sites(): 
                            site_exist_in_path = False 
                    if not site_exist_in_path: 
                        path_not_exist.append(path)     
            if not path_not_exist: 
                return False
            else: 
                print("Adding rings")
                add_path = random.choice(path_not_exist)
                for site in add_path:
                    ysite = site.pos[1]
                    if max([abs(ymin-ysite),abs(ymax-ysite)]) > self.ly:
                        print("Add site exceeds box contrains")
                        return False 
                for site in add_path:
                        self.syst[site] = self.onsite 
                        self.syst[self.mirror_sym_pairs[site]] = self.onsite
                self.syst[self.lat.neighbors()] = self.hop
                return True
        
        def _remove_ring(paths): 
            
            # Creating a temp copy in case something goes wrong during the 
            # ring removal processs
            temp_syst = copy.deepcopy(self.syst) 
            path_exist = []
            for path in paths: 
                print(len(path))
                if len(path) <7: 
                    site_exist_in_path = True 
                    for site in path: 
                        if site not in self.syst.sites(): 
                            site_exist_in_path = False 
                    if site_exist_in_path:
                        path_exist.append(path)      
            if not path_exist: 
                return False 
            elif len(path_exist) > 1: 
                print("Removing rings")        
                remove_path = random.choice(path_exist)
                symmetric_duplicates = []
                for site in remove_path: 
                    if self.mirror_sym_pairs[site] in remove_path: 
                        symmetric_duplicates.append(self.mirror_sym_pairs[site])
                for site in symmetric_duplicates: 
                    remove_path.remove(site)
                for site in remove_path: 
                    print(site.pos)
                    del self.syst[site]  
                    del self.syst[self.mirror_sym_pairs[site]]

                try: 
                    self.syst.eradicate_dangling()
                except: 
                    self.syst = copy.deepcopy(temp_syst) 
                    del temp_syst
                    return False 
                try:
                    if self._is_continous():
                        return True
                except: 
                    self.random_mirror_symmetric() 
                    return True 
                else: 
                    self.syst = copy.deepcopy(temp_syst) 
                    del temp_syst
                    return False 
            else: 
                return False

        moved = False 
       
        rand = random.uniform(0,1)
        if rand < 0.5: 
            move = _add_ring 
        else: 
            move = _remove_ring
        ntrail = 0
        while not moved:
            ntrail += 1
            [s,t] = random.choice(possible_head_tails)
            paths = nx.all_simple_paths(self.full_double_graph,s,t,cutoff=5)          
            moved = move(paths)
            if ntrail > 100:
                self.random_mirror_symmetric()
                if move == _add_ring:
                    move = _remove_ring
                else: 
                    move = _add_ring
            #if moved == False:
                #print(s.pos,t.pos)
                #print(self.full_double_graph.nodes[s],self.full_double_graph.nodes[t])
        
         
    def random_mirror_symmetric(self,symmetry=['mirror'],Ncentral=7): 
               
        min_width = get_width(Ncentral,self.lat)
        rect_L_plus_W = self.lx/4.0 + self.ly
        self.syst = None 
        
        def _get_random_2pts(self,L,w): 
            """ Used internally for randomly choosing 2pts so 
            the width of the structure at the boundary is atleast w
            """
            pt1 = random.uniform(0,self.ly)
            if 0 <= pt1 <= w: 
                pt2 = random.uniform(pt1+w,L)
            elif L - w <= pt1 <=L: 
                pt2 = random.uniform(0,pt1-w)
            else: 
                prob1 = (pt1-w)/(L-2*w)
                if random.random() <= prob1:
                    pt2 = random.uniform(0,pt1-w)
                else: 
                    pt2 = random.uniform(pt1+w,L)
            return pt1,pt2
        
        def _shape_from_lines(pos,offset,lineCoeff1,lineCoeff2):
                x,y = pos
                if 0<=x<self.lx/2:
                    val1 = np.polyval(lineCoeff1,abs(x))
                    val2 = np.polyval(lineCoeff2,abs(x))
                else:
                    val1 = np.polyval(lineCoeff1,abs(self.lx-x))
                    val2 = np.polyval(lineCoeff2,abs(self.lx-x))
                
                if 0.0<=y+offset<self.ly and  0<=x<self.lx:    
                    if val1<=y<=val2: 
                        return True 
                    else: 
                        return False
        
        ypt11,ypt12 = self._get_random_2pts(self.ly,min_width)
        ypt11,ypt12 = min(ypt11,ypt12),max(ypt11,ypt12)
        PTS1 = [[0,ypt11],[0,ypt12]]

                
        ypt21,ypt22 = self._get_random_2pts(rect_L_plus_W,min_width)
        ypt21,ypt22 = min(ypt21,ypt22),max(ypt21,ypt22)
        PTS2=[]
        for pt in [ypt21,ypt22]: 
            if pt > self.ly:
                PTS2.append([self.lx/2.0+self.ly-pt,self.ly])
            else: 
                PTS2.append([self.lx/2.0,pt])
            
        PTS2 = sorted(PTS2,key=lambda x: x[0],reverse = True)
        if abs(PTS2[0][0] - PTS2[1][0]) < 1.e-4: 
            PTS2 = sorted(PTS2,key=lambda x: x[1])
                
        lineCoeff1 =  np.polyfit([PTS1[0][0],PTS2[0][0]],
                                     [PTS1[0][1],PTS2[0][1]],1)
        
        lineCoeff2 =  np.polyfit([PTS1[1][0],PTS2[1][0]],
                                     [PTS1[1][1],PTS2[1][1]],1)
            #offset = (lineCoeff1[1]+lineCoeff2[1])/2
        offset = lineCoeff1[1]
        lineCoeff1[1] -= offset 
        lineCoeff2[1] -= offset
        
        syst = kwant.Builder(kwant.TranslationalSymmetry([self.lx,0]))
        syst[self.lat.shape((lambda pos: _shape_from_lines(pos,offset,
                                              lineCoeff1=lineCoeff1, 
                                              lineCoeff2=lineCoeff2)),(0,0))]=self.onsite 
        syst[self.lat.neighbors()]=self.hop
        syst.eradicate_dangling()
        self.syst = syst
    
    def random_inversion_symmetric(self,Ncentral=7): 
        min_width = get_width(Ncentral,self.lat)
        rect_L_plus_W = self.lx/4.0 + self.ly
        self.syst = None 
        
        def _get_randsym_2pts(L,w): 
            """ Used internally for randomly choosing 2pts so 
            the width of the structure at the boundary is atleast w
            """
            pt1 = random.uniform(0,L)
            pt2 = -1*random.uniform(0,L)
            if 0 <= abs(pt1)+abs(pt2) <= w:
                short_by = w - abs(pt1) - abs(pt2)
                pt1 += short_by/2.0 
                pt2 -= short_by/2.0
            return pt1,pt2
    
        def _shape_from_lines(pos,lineCoeff1,lineCoeff2):
                x,y = pos
                if -self.ly<=y<self.ly and  0<=x<self.lx: 
                    if 0<=x<self.lx/2:
                        val1 = np.polyval(lineCoeff1,abs(x))
                        val2 = np.polyval(lineCoeff2,abs(x))  
                        if val1<=y<=val2: 
                            return True 
                        else:
                            return False
                    else:
                        val1 = -1*np.polyval(lineCoeff1,abs(self.lx-x))
                        val2 = -1*np.polyval(lineCoeff2,abs(self.lx-x))
                        if -self.ly<=y<self.ly and  0<=x<self.lx:    
                            if val2<=y<=val1: 
                                return True
                            else: 
                                return False
                else: 
                    return False
                
        ypt11 = -1*random.uniform(0,self.ly)
        ypt12 = -1*ypt11
        if ypt12-ypt11 < min_width: 
            short_by = min_width - abs(ypt11)-abs(ypt12)
            ypt11 += -1*abs(short_by)
        #ypt11,ypt12 = min(ypt11,ypt12),max(ypt11,ypt12)
        PTS1 = [[0,ypt11],[0,ypt12]]
        ypt21,ypt22 = _get_randsym_2pts(rect_L_plus_W,min_width)
        ypt21,ypt22 = min(ypt21,ypt22),max(ypt21,ypt22)
        PTS2=[]
        for pt in [ypt21,ypt22]: 
            if pt > self.ly:
                PTS2.append([self.lx/2.0+self.ly-pt,self.ly])
            else: 
                PTS2.append([self.lx/2.0,pt])
            
        PTS2 = sorted(PTS2,key=lambda x: x[0],reverse = True)
        if abs(PTS2[0][0] - PTS2[1][0]) < 1.e-4: 
            PTS2 = sorted(PTS2,key=lambda x: x[1])
                
        lineCoeff1 =  np.polyfit([PTS1[0][0],PTS2[0][0]],
                                     [PTS1[0][1],PTS2[0][1]],1)
        
        lineCoeff2 =  np.polyfit([PTS1[1][0],PTS2[1][0]],
                                     [PTS1[1][1],PTS2[1][1]],1)
            #offset = (lineCoeff1[1]+lineCoeff2[1])/2
        #offset = lineCoeff1[1]
        #lineCoeff1[1] -= offset 
        #lineCoeff2[1] -= offset
        
        syst = kwant.Builder(kwant.TranslationalSymmetry([self.lx,0]))
        syst[self.lat.shape((lambda pos: _shape_from_lines(pos,
                                              lineCoeff1=lineCoeff1, 
                                              lineCoeff2=lineCoeff2)),(0,0))]=self.onsite 
        syst[self.lat.neighbors()]=self.hop
        syst.eradicate_dangling()
        self.syst = syst   
    
    def get_syst(self): 
        """
        Returns the kwant.Builder instance corresponding to the current state 
        of the random structure generator
        """ 
        return self.syst
    
    def set_syst(self,syst): 
        """
        Sets the current state of the random structure generator to the system 
        given 
        
        Parameters 
        ---------- 
        
        syst: kwant.Builder instance
        
        Returns 
        -------
            None 
        """
        self.syst=syst
    
    def syst2poscar(self,filename='POSCAR'):
        """
        Saves the current state of the random structure generator as a POSCAR file
        
        Parameters
        ----------
        
        filename: str
                Path to the POSCAR file
        """
        l = Lattice.from_lengths_and_angles([self.lx,self.ly,5.0],[90,90,90])
        pos = self._get_site_pos()
        nsites = len(pos)
        z = np.zeros(nsites)
        pos = np.column_stack((pos,z)) 
        structure = Structure(l,['C']*nsites,pos,coords_are_cartesian=True)
        structure.to(fmt='poscar',filename=filename)

    def poscar2syst(self,POSCAR): 
        """
        Reads a POSCAR file and generates a kwant system with the same geometery. 
        The onsite and hopping are those of the generator. lx, ly and lattice 
        of the generator must match that of the poscar files. Use 
        set_cell_attributes() to set lx,ly and lattice.  
        
        Parameters
        ---------- 
        
        POSCAR: str 
                path to POSCAR file 
        
        Returns 
        -------
            kwant.Builder() instance
        """
        struct = Structure.from_file(POSCAR)
        poscar_pos = np.array([[item.coords[0],item.coords[1]] for item in struct])
        min_x = np.min(poscar_pos[:,0])
        edge_sites = poscar_pos[poscar_pos[:,0]==min_x]
        min_y = np.min(edge_sites[:,1])
        poscar_pos[:,1] -= min_y
        syst = kwant.Builder(kwant.TranslationalSymmetry([self.lx,0]))
        def check_sites(pos):
            x,y = pos 
            for test_site in poscar_pos: 
                diff_x = abs(test_site[0]-x)
                diff_y = abs(test_site[1]-y)
                if diff_x < 1.0e-3 and diff_y < 1.0e-3 :
                    return True 
            return False
            
        syst[self.lat.shape(check_sites,(0,0))]=self.onsite
        syst[self.lat.neighbors()]=self.hop 
        self.syst = syst
        return self.syst 
    
        
    def translate_cell(self,t=None):  
        """Translates the x-coordinate of the sites by "t". 
        Can also be thought of as shifting/sliding the unit cell boundary over a fixed lattice
        
        Parameters
        ----------
        t = float 
            offset, default = lx/2.0 
        
        Returns 
        ------- 
            kwant.Builder() instance
        """
        if t is None: 
            t = self.lx/2.0 
        pos = np.array(self._get_site_pos())
        pos[:,0] += t
        for p in pos: 
            if p[0] >= self.lx: 
                p[0] -= self.lx
            elif p[0] <=0: 
                p[0] += self.lx 
        min_x = np.min(pos[:,0])
        edge_sites = pos[pos[:,0]==min_x]
        min_y = np.min(edge_sites[:,1])
        pos[:,1] -= min_y
        np.savetxt('pos.dat',pos)
        syst = kwant.Builder(kwant.TranslationalSymmetry([self.lx,0]))
        
        def _shape(site): 
            x,y = site
            for test_site in pos: 
                diff_x = abs(test_site[0]-x)
                diff_y = abs(test_site[1]-y)
                if diff_x < 1.0e-3 and diff_y < 1.0e-3 :
                    return True 
            return False

        syst[self.lat.shape(_shape,(0,0))] = self.onsite
        syst[self.lat.neighbors()]=self.hop
        self.syst = syst 
        return self.syst 
        
    def get_pol(self,wcc=False):
        """
        Returns the polarization of the geometry corresponding to the curent 
        state of the random structure generator
        
        Returns 
        -------
        
        pol: float 
            Polarization. Quantized at 0 or 0.5 if the system has spatial 
            symmetries
        
        """
        a = self.lx
        b = self.ly
        #for site in list(self.syst.sites(): 
        act_pos = np.array([site.pos for site in list(self.syst.sites())
                            if 0 <=site.pos[0]<=a])
        finalized_syst = self.syst.finalized()
        red_pos = np.zeros(np.shape(act_pos))
        red_pos[:,0] = act_pos[:,0]/a
        red_pos[:,1] = act_pos[:,1]/b
        return calc_pol(finalized_syst,red_pos,wcc=wcc);
    
    def plot_syst(self): 
        kwant.plot(self.syst)
        
    def _get_bands(self,momenta=(-np.pi,np.pi)): 
        bands = kwant.physics.Bands(self.syst.finalized())
        momenta = np.linspace(momenta[0],momenta[1],101)
        eigs = [bands(k,return_eigenvectors=True) for k in momenta]
        return eigs 
    
    def _get_density(self): 
        eigs = self._get_bands() 
        eigvals,eigvecs = zip(*eigs)
        eigvecs = np.array(eigvecs)
        nbands = len(eigvals)
        n_occupied_bands = int(nbands/2)
        k_sum_states = np.sum(eigvecs[:,:,:],axis=0) # numbands = numbasis 
        sum_states = np.sum(k_sum_states[:,:n_occupied_bands],axis=1)
        #sum_states = np.sum(k_sum_states[:,:],axis=1)
        density = np.absolute(sum_states)
        print(density)
        return density
    
    def _get_density2(self): 
        eigs = self._get_bands() 
        eigvals,eigvecs = zip(*eigs)
        eigvecs = np.array(eigvecs)
        nbands = len(eigvals)
        n_occupied_bands = int(nbands/2)
        nbasis = np.shape(eigvecs)[-1]
        coefficients = np.zeros([nbasis])
        for k in range(len(eigvecs)): 
            states = np.absolute(eigvecs[k,:,:n_occupied_bands])
            sum_states = np.zeros([nbasis])
            for i in range(n_occupied_bands):
                sum_states += states[:,i]/np.sum(states[:,i])
            coefficients += sum_states/(len(eigvecs)*n_occupied_bands)
        #k_sum_states = np.sum(eigvecs[:,:,:],axis=0) # numbands = numbasis 
        #sum_states = np.sum(k_sum_states[:,:n_occupied_bands],axis=1)
        #sum_states = np.sum(k_sum_states[:,:],axis=1)
        #density = np.absolute(sum_states)
        print(coefficients)
        print(sum(coefficients))
        return coefficients     
        
    
    def _get_density3(self): 
        eigs = self._get_bands() 
        eigvals,eigvecs = zip(*eigs)
        eigvecs = np.array(eigvecs)
        nbands = len(eigvals)
        n_occupied_bands = int(nbands/2)
        nbasis = np.shape(eigvecs)[-1]
        rho = kwant.operator.Density(self.syst.finalized())
        s = np.zeros([nbasis])
        for i in range(len(eigvecs)): 
            for j in range(n_occupied_bands): 
               s = rho(eigvecs[i,:,j])
        return s
     
    def _1D_to_finite(self): 
        pos_lattice = np.array(self._get_site_pos())
        syst = kwant.Builder()
        def check_sites(pos):
            x,y = pos 
            for test_site in pos_lattice: 
                diff_x = abs(test_site[0]-x)
                diff_y = abs(test_site[1]-y)
                if diff_x < 1.0e-3 and diff_y < 1.0e-3 :
                    return True 
            return False
            
        syst[self.lat.shape(check_sites,(0,0))]=self.onsite
        syst[self.lat.neighbors()]=self.hop 
        return syst.finalized() 
    
    def _plot_density(self): 
        density = self._get_density3()
        fig = kwant.plotter.density(self._1D_to_finite(),density, relwidth=0.08,cmap='jet',background='white')#,oversampling=12);
        #return fig

        
    def get_adjacency(self): 
        """ 
        Returns the Adjacency Matrix.
        The adjacency matrix is built from the hopping pairs of the kwant system 
        so it is symmetric. 
   
        
        Returns 
        ------- 
        adjMat: numpy array 
                Symmetric Adjacency matrix 

        """ 
        nmax_sites = len(self.full_syst_pos)
        adjMat = np.zeros((nmax_sites,nmax_sites))   
        
        pos = np.array(self._get_site_pos())
        ymin = np.min(pos[:,1]) 
        offset = copysign(1,ymin)*ceil(abs(ymin)/self.lat.prim_vecs[1][1])*self.lat.prim_vecs[1][1]
                
        def _get_index(site):
            x,y = site 
            x -= floor(x/self.lx)*self.lx
            y -= offset
            for i,item in enumerate(self.full_syst_pos):
                diff = abs(item[0]-x) + abs(item[1]-y)
                if diff < 1.e-2:
                    return i
            return None 
        
        def _is_inside_cell(site1,site2): 
            x1,y1 = site1.pos 
            x2,y2 = site2.pos 
            if 0<=x1<=self.lx and 0<=x2<=self.lx: 
                return True 
            else: 
                return False

        for hopping, value in self.syst.hopping_value_pairs(): 
            site1,site2=hopping
            index_i = _get_index(site1.pos)
            index_j = _get_index(site2.pos)
            if _is_inside_cell(site1,site2):
                adjMat[index_i,index_j]=1
            else: 
                adjMat[index_i,index_j]=-1      
        return adjMat.astype('int8')

    def construct_graph(self,draw=False): 
        """
        Construct a directed graph from the kwant system. The direction of the 
        graph is from left to right (increasing x). 
        
        Parameters 
        ----------
        draw: boolean 
            If true, outputs visualization of the graph 
        
        Returns 
        -------
        G : nx.DiGraph instance 
            Graph of the kwant system
        """
        
        G = nx.DiGraph()
        adjMat = self.get_adjacency() 
        I,J = np.nonzero(adjMat)
        for i,j in zip(I,J):
            if adjMat[i,j] > 0: 
                hop=1
            else: 
                hop=-1
            x1,y1 = self.full_syst_pos[i]
            x2,y2 = self.full_syst_pos[j] 
            G.add_node(i,x=x1,y=y1)
            G.add_node(j,x=x2,y=y2)
            if x1 < x2: 
                G.add_edge(i,j,hop=hop)
            else: 
                G.add_edge(j,i,hop=hop)
        self.graph = G
        if draw: 
            self.draw_lattice_graph()
        return G 

    def draw_lattice_graph(self): 
        """
        Draws the kwant system accroding to the position of the sites
        
        Returns 
        -------
        None 
        """
        pos = {}        
        for node in self.graph.nodes(data=True):
            pos[node[0]] = [node[-1]['x'],node[-1]['y']]
        edge_color=[]
        for edge in self.graph.edges(data=True): 
            if edge[-1]['hop'] > 0: 
                edge_color.append('black')
            elif edge[-1]['hop'] < 0: 
                edge_color.append('red')
        nx.draw_networkx(self.graph,pos=pos,edge_color=edge_color)
    
    def _find_edge_connections(self):
        edge_connections = []
        for edge in self.graph.edges(data=True):
            if edge[-1]['hop'] < 0: 
                edge_connections.append([edge[0],edge[1]])
        return edge_connections
    
    def find_paths_across_cell(self):
        """
        Find the shortest simple paths from the left end of the unit cell to 
        the right end for each bond crossing the unit cell 
        
        Returns 
        -------
        paths: dict 
                key: index of head, index of tail 
                values: paths 
        """
        if self.graph is None: 
            G=self.construct_graph()
        edge_connections = self._find_edge_connections()
        #print(edge_connections)
        paths = {}
        for connection in edge_connections:
            short_paths = []
            for p in nx.shortest_simple_paths(self.graph,connection[0],connection[1]):
                short_paths.append(p) 
            paths[connection[0],connection[1]]=short_paths
        return paths
        
        
        
    
                        
            
                
        
    
        
