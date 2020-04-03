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
from math import floor

class EncodeLattice():
    def __init__(self,lat='Armchair',Nlx=10,Nly=20):
        self.lat = lattices[lat]
        self.lx = self.lat.prim_vecs[0][0]
        self.ly = self.lat.prim_vecs[1][1]
        self.UnitCell = self.getUnitCell(self.lat)
        self.UnitCellPos = self._get_pos_list(self.UnitCell)
        self.encodeDict = self._get_encode_dict()
        self.nly = Nly         
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
        # Shifting the lattice if there are sites with y < 0 
        min_y = np.min(pos_list[:,1])
        if min_y < 0: 
            shiftVector = floor(min_y/self.ly)*self.ly
            pos_list[:,1] -= shiftVector
        nlx = int(round(np.max(pos_list[:,0])/self.lx))
        nly = self.nly  
        systCode=[]

        for ny in range(nly): 
            rowCode=[]
            for nx in range(nlx): 
                x0 = self.lx*nx 
                y0 = self.ly*ny 
                shiftedLat = pos_list - [x0,y0]
                mask = (0<=shiftedLat[:,0]) & (shiftedLat[:,0]<self.lx) & \
                        (0<=shiftedLat[:,1]) & (shiftedLat[:,1]<self.ly-0.05) # 0.05 is to ensure the condition holds true 
                                                                              # even when there are numerical errors 
                CompareUnit = shiftedLat[mask,:]
                codeFound=False
                for code,lat in self.encodeDict.items():
                    latArray = np.array(lat)
                    if latArray.shape == CompareUnit.shape:
                        boolArray = []
                        for site in latArray: 
                            boolArray.append(self._site_is_in(CompareUnit,site))
                        if np.all(np.array(boolArray)): 
                            tileCode = code 
                            codeFound = True 
                if codeFound==False: 
                    tileCode=0
                rowCode.append(tileCode)
            systCode.append(rowCode)
        return np.flip(np.array(systCode),axis=0)
                            
                    
                    
                
                

        

