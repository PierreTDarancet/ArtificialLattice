#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 16:27:19 2019

@author: ssrinivasan
"""
#from helper import *
from .utilities import lattices,get_width
from .dump import dump 
from pymatgen import Lattice, Structure
from math import floor,ceil,copysign
from operator import itemgetter
from io import StringIO
import kwant,z2pack,cmath,random,operator,copy,itertools 
import numpy as np 
import networkx as nx 
import matplotlib.pyplot as plt 
#import networkx.algorithms.isomorphism as iso

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
                              pos_tol=1e-2,iterator=range(7,1001,2));
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
    
    def __init__(self,lat='graphene',nlx=5,
                 nly= 5,onsite=0,hop=-1,syst=None,dumpfile=None): 
        self.lat = lattices[lat]
        self.lx  = nlx*self.lat.prim_vecs[0][0]
        self.ly  = nly*self.lat.prim_vecs[1][1]
        self.onsite = onsite
        self.hop = hop
        self.syst = syst
        self.full_double_syst_pos = None
        self.index_dict = None 
        self.full_graph = None
        self.full_double_graph = None 
        self.full_syst = None
        self.full_double_syst = None 
        self.sym_pairs = None 
        self.full_syst_lat = None 
        self._set_full_syst_attributes()
        self.graph = None 
        
        
        if dumpfile:
            self.dumpfile = dumpfile 
            global d 
            d = dump(filename=dumpfile)
        else: 
            self.dumpfile = False        
        
    def _get_mirror_symmetric_pairs(self): 
        mirror_plane = self.lx/2.0
        mirror_sym_pairs = {}
        for site1 in list(self.full_double_syst.sites()): 
            x1,y1 = site1.pos
            for site2 in list(self.full_double_syst.sites()): 
                if site1 is not site2: 
                    x2,y2 = site2.pos
                    x_check = np.isclose(abs(x1-mirror_plane),abs(x2-mirror_plane))  
                    y_check = np.isclose(y1,y2)
                    if x_check and y_check:
                        mirror_sym_pairs[site1] = site2
        return mirror_sym_pairs
    
    def _get_inversion_symmetric_pairs(self): 
        inversion_plane = self.lx/2.0
        inv_sym_pairs = {}
        for site1 in list(self.full_double_syst.sites()): 
            x1,y1 = site1.pos
            for site2 in list(self.full_double_syst.sites()): 
                if site1 is not site2: 
                    x2,y2 = site2.pos
                    x_check = np.isclose((x1-inversion_plane),-1*(x2-inversion_plane))  
                    y_check = np.isclose(y1,-1*y2)
                    if x_check and y_check:
                        inv_sym_pairs[site1] = site2
        return inv_sym_pairs
        
                                                
    def _get_site_pos(self,syst=None):
        if syst is None:
            pos = [site.pos for site in list(self.syst.sites())]
            return pos
        else: 
            pos = [site.pos for site in list(syst.sites())]
            return pos
        
    def _set_full_syst_attributes(self): 
        self.full_syst = self._make_full_syst(uly=0)
        self.full_graph = self._construct_full_graph(syst=self.full_syst)
        self.full_double_syst = self._make_full_syst(uly=-1*self.ly) 
        pos = np.array(sorted(self._get_site_pos(syst=self.full_double_syst),
                     key=operator.itemgetter(0,1)))
        index_dict = {}
        for i,site in enumerate(pos): 
            index_dict[tuple(site)] =i 
        self.index_dict = index_dict
        self.full_double_syst_pos = pos
        self.full_double_graph = self._construct_full_graph(syst=self.full_double_syst)
        self._fill_all_sites()
        #self.graph = self.construct_graph()
        self.sym_pairs = {'mirror':self._get_mirror_symmetric_pairs(),
                          'inversion':self._get_inversion_symmetric_pairs()}
        return None 
    
    def _make_full_syst(self,uly=0,trans_sym_direction='x'):
        def _template():
            syst = kwant.Builder()
            syst[self.lat.shape((lambda pos: 0<=pos[0]<self.lx  \
                and 0<=pos[1]<self.ly),(0,0))] = self.onsite 
            syst[self.lat.neighbors()] = self.hop 
            return syst
        pos = self._get_site_pos(syst=_template())
        if trans_sym_direction=='x':
            a = self.lx
            b = self.ly
            trans_vec=[a,0]
        elif trans_sym_direction=='y':
            a = self.lx+0.1
            b = self.ly
            trans_vec=[0,b]
        else:
            raise #"Translation Symmetry direction should be 'x' or 'y'"
        lattice_1D = kwant.lattice.general([[a,0],[0,b]],pos)
        self.full_syst_lat = lattice_1D
        syst_1D = kwant.Builder(kwant.TranslationalSymmetry(trans_vec))
        if trans_sym_direction=='x':
            syst_1D[lattice_1D.shape((lambda pos: uly<= pos[1] < self.ly),(0,0))]=self.onsite
        if trans_sym_direction=='y':
            syst_1D[lattice_1D.shape((lambda pos: 0<= pos[0] < a),(0,0))]=self.onsite
        syst_1D[lattice_1D.neighbors()] = self.hop
        return syst_1D
    
    def set_cell_attributes(self,nlx=5,nly=5,lat='Armchair'): 
        self.lat = lattices[lat]
        self.lx = nlx*self.lat.prim_vecs[0][0] 
        self.ly = nly*self.lat.prim_vecs[1][1] 
        self._set_full_syst_attributes() 
    
    def _fill_all_sites(self): 
        full_syst=self._make_full_syst() 
        self.syst = full_syst


    def _finite_to_1D(syst,lat_vec,trans_sym_direction='x'): 
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

        sites = list(syst.sites())
        pos = [site.pos for site in sites]
        if trans_sym_direction=='x':
            a = lat_vec
            b = max(np.array(pos)[:,1])+2.0
            trans_vec=[a,0]
        elif trans_sym_direction=='y':
            a = max(np.array(pos)[:,0])+2.0
            b = lat_vec
            trans_vec=[0,b]
        else:
            raise #"Translation Symmetry direction should be 'x' or 'y'"
        lattice_1D = kwant.lattice.general([[a,0],[0,b]],pos)
        syst_1D = kwant.Builder(kwant.TranslationalSymmetry(trans_vec))
        if trans_sym_direction=='x':
            syst_1D[lattice_1D.shape((lambda pos: 0< pos[1] <= b),(0,0))]=0
        if trans_sym_direction=='y':
            syst_1D[lattice_1D.shape((lambda pos: 0< pos[0] <= a),(0,0))]=0
        syst_1D[lattice_1D.neighbors()] = -1
        return syst_1D

        
    def _construct_full_graph(self,draw=False,syst=None):
        if not syst: 
            syst = self.syst
        G = nx.Graph()
        for hopping,value in syst.hopping_value_pairs():
            u,v = hopping
            G.add_node(u,x=u.pos[0],y=u.pos[1]) 
            G.add_node(v,x=v.pos[0],y=v.pos[1]) 
            G.add_edge(u,v,hop=value)
        return G 
        
    def _is_continous(self): 
        G = self.construct_graph() 
        G_undir = G.to_undirected()
        graph_xlist = [node for node in G.nodes('x')]
        graph_xlist.sort(key=itemgetter(1))
        head = graph_xlist[0][0]
        tailx= graph_xlist[-1][1] 
        head_neigh = G.neighbors(head)
        is_continous = False
        for node in head_neigh: 
            if abs(G.nodes[node]['x']-tailx) < 1.e-3: 
                is_continous = True
        remove_edges = []
        for edge in G_undir.edges(data=True): 
            if edge[-1]['hop'] < 0: 
                remove_edges.append([edge[0],edge[1]])
        for edge in remove_edges: 
            G_undir.remove_edge(edge[0],edge[1])
        is_cluster = nx.is_connected(G_undir)
        return is_cluster and is_continous
    
    
    def _get_possible_head_tails(self): 
        edge_sites = []
        for site in self.syst.sites(): 
            if self.syst.degree(site)<3:
                edge_sites.append(site.pos)
        node_list = []
        for site in edge_sites:
            x,y = site
            for node in self.full_double_graph.nodes(data=True): 
                if np.isclose(node[-1]['x'],x) and np.isclose(node[-1]['y'],y): 
                    node_list.append(node[0])
                    
        possible_head_tails = []
        for (s,t) in itertools.combinations(node_list,2):
            if nx.shortest_path_length(self.full_double_graph.to_undirected(),
                                       s,t) <7:
                possible_head_tails.append([s,t])
        return possible_head_tails
    
    def swap_move(self,sym='mirror'):
        if sym not in ['mirror','inversion']: 
            raise('sym argument needs to be either "mirror" or "inversion"')
        if sym =='mirror':
            lines_move = self.random_mirror_symmetric
        else: 
            lines_move = self.random_inversion_symmetric
        sym_pairs = self.sym_pairs[sym]  
        possible_head_tails = self._get_possible_head_tails()
        
        def _add_ring(paths):
            width = 0
            temp_syst = copy.deepcopy(self.syst) 
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
                    w = max([abs(ymin-ysite),abs(ymax-ysite)])
                    width = max(width,w)
                    if width > self.ly: 
                        self.syst = copy.deepcopy(temp_syst) 
                        del temp_syst
                        return False
                    else: 
                        self.syst[site] = self.onsite 
                        self.syst[sym_pairs[site]] = self.onsite
                self.syst[self.full_syst_lat.neighbors()] = self.hop
                return True
        
        def _remove_ring(paths): 
            print("Removing rings") 
            temp_syst = copy.deepcopy(self.syst) 
            path_exist = []
            for path in paths: 
                if len(path) <7: 
                    site_exist_in_path = True 
                    for site in path: 
                        print(site.pos)
                        if site not in self.syst.sites(): 
                            site_exist_in_path = False 
                    if site_exist_in_path:
                        path_exist.append(path)   
            print(path_exist)
            if not path_exist: 
                return False 
            else:
                print("Removing rings")        
                remove_path = random.choice(path_exist)
                symmetric_duplicates = []
                for site in remove_path: 
                    if sym_pairs[site] in remove_path: 
                        symmetric_duplicates.append(sym_pairs[site])
                for site in symmetric_duplicates: 
                    remove_path.remove(site)
                for site in remove_path: 
                    #print(site.pos)
                    del self.syst[site]  
                    del self.syst[sym_pairs[site]]
                    #kwant.plot(self.syst)

                try: 
                    self.syst.eradicate_dangling()
                except: 
                    self.syst = copy.deepcopy(temp_syst) 
                    del temp_syst
                    print("Eradicate fail")
                    return False 

                if self._is_continous():
                    return True
                else: 
                    #kwant.plot(self.syst)
                    #self.draw_lattice_graph()
                    self.syst = copy.deepcopy(temp_syst) 
                    del temp_syst
                    print("Continue fail")
                    return False 


        moved = False 
        rand = random.uniform(0,1)
        if rand < 0.5: 
            move = _add_ring 
        else: 
            move = _remove_ring
        ntrail = 0
        kwant.plot(self.syst)
        while not moved:
            ntrail += 1
            [s,t] = random.choice(possible_head_tails)
            print(s.pos,t.pos)
            paths = nx.all_simple_paths(self.full_double_graph.to_undirected()
                                        ,s,t,cutoff=5)          
            moved = move(paths)
            if ntrail > 10:
                lines_move()
                possible_head_tails = self._get_possible_head_tails()
                if move == _add_ring:
                    move = _remove_ring
                else: 
                    move = _add_ring        
     
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
        
    def random_mirror_symmetric(self,symmetry=['mirror'],Ncentral=7): 
               
        min_width = get_width(Ncentral,self.lat)
        rect_L_plus_W = self.lx/4.0 + self.ly
        self.syst = None 
    
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
        
        syst = kwant.Builder(kwant.TranslationalSymmetry([self.lx,0]))
        syst[self.lat.shape((lambda pos: _shape_from_lines(pos,
                                              lineCoeff1=lineCoeff1, 
                                              lineCoeff2=lineCoeff2)),(0,0))]=self.onsite 
        syst[self.lat.neighbors()]=self.hop
        temp_syst = copy.deepcopy(syst)
        try:
            temp_syst.eradicate_dangling()
            syst = copy.deepcopy(temp_syst)
        except: 
            kwant.plot(syst)
            raise
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
        
    def syst2dump(self,frame=0): 
        if self.dumpfile:
            l = Lattice.from_lengths_and_angles([self.lx,self.ly,5.0],[90,90,90])
            lammpsBox=[[0.0,l.matrix[0,0]],[0.0,l.matrix[1,1]],
                       [0.0,l.matrix[2,2]]]
            timestep = frame 
            atoms_data={} 
            pos = self._get_site_pos()
            natoms = len(pos)
            z = np.zeros(natoms)
            id_col = np.arange(1,natoms+1)
            print("Writing frame {} to dump file".format(frame))
            atoms_data['attributes']=['id','x','y','z']
            data = np.column_stack((id_col,pos,z))
            atoms_data['data'] = data
            d.write_frame(timestep,lammpsBox,natoms,atoms_data)
        else: 
            raise('StructGen was not declared with a dump file name')
    
    def dump2syst(self,frame): 
        def _read_frame(frame):
            xbox,ybox =frame.split('\n',9)[-5:-3]
            lx,ly = float(xbox.split()[-1]),float(ybox.split()[-1])
            xy = np.loadtxt(StringIO(frame),skiprows=9,usecols=(1,2))
            return [lx,ly,xy]
        lx,ly,xy = _read_frame(frame)
        min_x = np.min(xy[:,0])
        edge_sites = xy[xy[:,0]==min_x]
        min_y = np.min(edge_sites[:,1])
        xy[:,1] -= min_y
        lat_y = 2*ly+2.0
        lat = kwant.lattice.general([[lx,0],[0,lat_y]],xy)
        syst = kwant.Builder(kwant.TranslationalSymmetry([lx,0]))
        syst[lat.shape((lambda pos: 0 <= pos[1] < lat_y),(0,0))]= self.onsite
        syst[lat.neighbors()] = self.hop
        self.syst = syst 
        return syst 
    
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
        nbasis = np.shape(eigvecs)[-1]
        nbands = nbasis
        n_occupied_bands = int(nbands/2)
        rho = kwant.operator.Density(self.syst.finalized())
        s = np.zeros([nbasis])
        for i in range(len(eigvecs)): 
            for j in range(n_occupied_bands): 
               s += rho(eigvecs[i,:,j])
        return s/(len(eigvecs)*n_occupied_bands)
     
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
        density = self._get_density()
        fig = kwant.plotter.map(self._1D_to_finite(),density,oversampling=50)
        
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
        nmax_sites = len(self.full_double_syst_pos)
        adjMat = np.zeros((nmax_sites,nmax_sites))   
        pos = np.array(self._get_site_pos())
        ymin,ymax = np.min(pos[:,1]),np.max(pos[:,1]) 
        clearance = self.ly - (ymax-ymin)
        #clearance = ymin
        if clearance <= 0: 
            offset = copysign(1,clearance)*(ceil(abs(clearance)/self.lat.prim_vecs[1][1])+1)*self.lat.prim_vecs[1][1]
        else: 
            offset=0
        print('offset={}'.format(offset))
        def _get_index(site):
            x,y = site 
            # Remap the sites outsite the box (images) to the corresponding
            # sites inside box 
            x -= floor(x/self.lx)*self.lx
            y += offset
            for i,item in enumerate(self.full_double_syst_pos):
                diff = abs(item[0]-x) + abs(item[1]-y)
                if diff < 1.e-2:
                    return i
            print('Matching node for {},{} not found in double_syst_graph'.format(x,y))
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
                if None not in [index_i,index_j]: adjMat[index_i,index_j]=1
                else:
                    #self.draw_lattice_graph()
                    self.plot_syst()
            else: 
                if None not in [index_i,index_j]: adjMat[index_i,index_j]=-1
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
            x1,y1 = self.full_double_syst_pos[i]
            x2,y2 = self.full_double_syst_pos[j] 
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

    def draw_lattice_graph(self,graph=None,figsize=None): 
        """
        Draws the kwant system accroding to the position of the sites
        
        Returns 
        -------
        None 
        """
        if not graph: 
            graph = self.construct_graph()
        pos = {}        
        for node in graph.nodes(data=True):
            pos[node[0]] = [node[-1]['x'],node[-1]['y']]
        edge_color=[]
        width=[]
        for edge in graph.edges(data=True): 
            if edge[-1]['hop'] > 0: 
                edge_color.append('black')
                width.append(1.0)
            elif edge[-1]['hop'] < 0: 
                edge_color.append('red')
                width.append(0.5)
        if figsize is None: 
            figsize = (8,6)
        nx.draw_networkx(graph,pos=pos,edge_color=edge_color,node_size=500,
                         alpha=0.7,width=1.4,with_labels=False,linewidths=2.0)
    
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
        paths = {}
        for connection in edge_connections:
            short_paths = []
            for p in nx.shortest_simple_paths(self.graph,connection[0],connection[1]):
                short_paths.append(p) 
            paths[connection[0],connection[1]]=short_paths
        return paths
        
    def terminate_edges(self,bulk_degree=3,delta_t=-0.2):
        sites = [ site for site in self.syst.sites()]
        edges = []
        for site in sites:
            if self.syst.degree(site) < bulk_degree: 
                edges.append(site)
                for neigh in self.syst.neighbors(site): 
                    self.syst[site,neigh] += delta_t/2.0 
        
                
        
        
        
    
                        
            
                
        
    
        
