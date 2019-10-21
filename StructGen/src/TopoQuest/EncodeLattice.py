#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 13:13:23 2019

@author: ssrinivasan
"""

from TopoQuest.utilities import lattices 
from itertools import combinations
import numpy as np 
import kwant 


class EncodeLattice():

    def __init__(self,lat='Armchair',Nlx=10,Nly=10):
        self.lat = lattices[lat]
        self.lx = self.lat.prim_vecs[0][0]
        self.ly = self.lat.prim_vecs[1][1]
        self.UnitCell = self.getUnitCell(self.lat)
        self.UnitCellPos = self._get_pos_list(self.UnitCell)
        self.center = [self.lx/2.0, self.ly/2.0]
        self.encodeDict = self._get_encode_dict()
        
    def _get_encode_dict(self): 
        encodeDict={}
        vacancy_list = range(len(self.UnitCellPos)+1)
        i=0
        for v in vacancy_list:
            for tile in combinations(self.UnitCellPos,v): 
                encodeDict[i]=tile 
                i += 1
        return encodeDict
        
        
    def _get_lat_center(lat): 
        prim_vecs = lat.prim_vecs
        xc = prim_vecs[0][0]/2.0 
        yc = prim_vecs[1][1]/2.0 
        return [xc,yc]
    
    def getUnitCell(self,lat): 
        """Returns a list of positions of the unit cell for a given lattice
            (Currently works only for orthogonal unit cell) 
        """ 
        
        # TODO: Works only for orthogonal unit cell. 
        # Find solutions for triclinic ones 
        prim_vecs = lat.prim_vecs
        a = prim_vecs[0][0] 
        b = prim_vecs[1][1]
        
        unitCell = kwant.Builder()
        unitCell[lat.shape((lambda pos: 0<=pos[0]<a and 0<=pos[1]<b),(0,0))] = 0 
        unitCell[lat.neighbors()]=-1 
        return unitCell
    
    
    
    def _get_pos_list(self,syst):
        pos_list = []
        for site in syst.sites():
            pos_list.append(site.pos)
        return pos_list 
    
    
    def _site_is_in(self,syst_pos,site): 
        return np.isclose(syst_pos-site,[0,0]).all(axis=1).any()
    
    def encode(self,syst):
        pos_list = np.array(self._get_pos_list(syst))
        nlx = int(round(np.max(pos_list[:,0])/self.lx))
       # nly = int(round(np.max(pos_list[:,1]) - np.min(pos_list[:,1])/self.ly) + 1)
        nly = 10  
        systCode=[]
        for ny in range(nly): 
            rowCode=[]
            for nx in range(nlx): 
                xc = self.center[0] +self.lx*nx 
                yc = self.center[1] +self.ly*ny 
                shiftedLat = pos_list - [xc,yc] 
                codeFound=False
                for code,lat in self.encodeDict.items():
                    if code!=0:
                        for site in np.array(lat): 
                            if not self._site_is_in(shiftedLat,site): 
                                break 
                            else: 
                                tileCode=code
                                codeFound=True 
                if codeFound==False: 
                    tileCode=0
                rowCode.append(tileCode)
            systCode.append(rowCode)
            
        return np.flip(np.array(systCode),axis=0)
                            
                    
                    
                
                

        

