#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 16:27:19 2019

@author: ssrinivasan
"""

#from helper import *
from helper import Armchair,get_width
from pymatgen import Lattice, Structure
from math import floor
import kwant,z2pack
import cmath,random,operator 
import numpy as np 


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


def calc_pol(syst,red_pos):   
    """ Returns the polarization of a unit cell built using kwant
    
    Parameters:
    -----------
    syst: finalized kwant builder 
    red_pos: reduced position of the sites in the unit cell
    
    Returns: 
    --------
    result.pol: Polarization from the result of the line calculation using Z2pack
    """
    
    Hk = get_Hk(syst,dim=2)
    z2_system = z2pack.hm.System(Hk,dim=2,pos=red_pos,
                                 convention=2)
    result = z2pack.line.run(system=z2_system, 
                              line=lambda t1: [t1,0],
                              pos_tol=1e-2,iterator=range(7,501,2));
    return result.pol


class Check_redundant(): 
    """ Class containing methods to check for redundancy between the randomly
    generated structures by StructGen()
    
    Attributes: 
    ----------- 
    seen_lat: set of set of sites of previously seen unit cells  
    """
    def __init__(self):
        self.seen_lat = set()
    
    def is_redundant(self,syst):
        """Returns True if the structure has been seen before by the generator
        
        Parameters:
        -----------
        syst: Predicted builder by the generator 
        
        Returns: 
        --------
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
        
CR = Check_redundant()     

class StructGen(): 
    """
    Class for random 1D structure generator 
    
    Attributes: 
    -----------
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
        self.full_syst_pos,self.index_dict = self._index_sites()
        
    def _get_site_pos(self,syst=None):
        if syst is None:
            pos = np.array([site.pos for site in list(self.syst.sites())]).tolist()
            return pos
        else: 
            pos = np.array([site.pos for site in list(syst.sites())]).tolist()
            return pos
        
    def _index_sites(self): 
        full_syst = kwant.Builder()
        full_syst[self.lat.shape(
                (lambda pos: 0<=pos[0]<=self.lx and -1*self.ly<=pos[1]<=self.ly),
                (0,0))] = self.onsite 
        full_syst[self.lat.neighbors()] = self.hop
        pos = self._get_site_pos(syst=full_syst)
        pos.sort(key=operator.itemgetter(0,1))
        pos = np.array(pos)
        index_dict = {}
        for i,site in enumerate(pos): 
            index_dict[tuple(site)] =i 
        del full_syst 
        return pos,index_dict
        
        
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
        ntrial=0 
        while CR.is_redundant(self.syst):
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
            ntrial +=1 
            if ntrial > 10000: 
                return False 
        
        return True 
        
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
        
        Parameters: 
        ---------- 
        
        syst: kwant.Builder instance
        
        Returns: 
        ---------
            None 
        """
        self.syst=syst
    
    def syst2poscar(self,filename='POSCAR'):
        """
        Saves the current state of the random structure generator as a POSCAR file
        
        Parameters:
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
        The onsite and hopping are those of the generator
        
        Parameters:
        ---------- 
        
        POSCAR: str 
                path to POSCAR file 
        
        Returns: 
        ---------
            kwant.Builder() instance
        """
        struct = Structure.from_file(POSCAR)
        poscar_pos = np.array([[item.coords[0],item.coords[1]] for item in struct])
        min_x = np.min(poscar_pos[:,0])
        edge_sites = poscar_pos[poscar_pos[:,0]==min_x]
        min_y = np.min(edge_sites[:,1])
        poscar_pos[:,1] -= min_y
        self.lx = struct.lattice.a 
        self.ly = struct.lattice.b 
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
        syst.eradicate_dangling() 
        self.syst = syst 
        return syst 
        
    def translate_cell(self,t=None):  
        """Translates the x-coordinate of the sites by "t". 
        Can also be thought of as shifting/sliding the unit cell boundary over a fixed lattice
        
        Parameters:
        -----------
        t = float 
            offset, default = lx/2.0 
        
        Returns: 
        ---------- 
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
        self.lat = kwant.lattice.Polyatomic([[self.lx,0],[0,self.ly]],pos)
        syst = kwant.Builder(kwant.TranslationalSymmetry([self.lx,0]))
        syst[self.lat.shape((lambda pos: pos[1]>=0 and pos[1]<=self.ly),
                            (0,0))] = self.onsite
        syst[self.lat.neighbors()]=self.hop
        self.syst = syst 
        return syst 
        
    def get_pol(self):
        """
        Returns the polarization of the geometry corresponding to the curent 
        state of the random structure generator
        
        Returns: 
        -------
        
        pol: float 
            Polarization. Quantized at 0 or 0.5 if the system has spatial 
            symmetries
        
        """
        a = self.lx
        b = self.ly
        #for site in list(self.syst.sites(): 
        act_pos = np.array([site.pos for site in list(self.syst.sites()) \
                            if 0 <=site.pos[0]<=a])
        finalized_syst = self.syst.finalized()
        red_pos = np.zeros(np.shape(act_pos))
        red_pos[:,0] = act_pos[:,0]/a
        red_pos[:,1] = act_pos[:,1]/b
        return calc_pol(finalized_syst,red_pos);
    
    def plot_syst(self): 
        kwant.plot(self.syst)
    

    def get_adjacency(self): 
        """ 
        Returns the Adjacency Matrix.
        The adjacency matrix is built from the hopping pairs of the kwant system 
        so it is symmetric. 
   
        
        Returns: 
        ------- 
        adjMat: numpy array 
                Symmetric Adjacency matrix 

        """ 
        nmax_sites = len(self.full_syst_pos)
        adjMat = np.zeros((nmax_sites,nmax_sites))   
        
        def _get_index(site):
            x,y = site 
            x -= floor(x/self.lx)*self.lx
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
        return adjMat
                
            
                
        
    
        
