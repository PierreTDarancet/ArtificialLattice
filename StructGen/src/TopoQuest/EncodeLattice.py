#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 13:13:23 2019

@author: ssrinivasan
"""

from TopoQuest.utilities import lattices 
import kwant 


class EncodeLattice():

    def __init__(self,lat='Armchair',Nlx=10,Nly=10):
        self.lat = lattices[lat]
        self.lx = self.lat.prim_vecs[0][0]
        self.ly = self.lat.prim_vecs[1][1]
        self.UnitCellPosList = self.getUnitCellPos(self.lat)
        self.center = [self.lx/2.0, self.ly/2.0]
        
    def _get_lat_center(lat): 
        prim_vecs = lat.prim_vecs
        xc = prim_vecs[0][0]/2.0 
        yc = prim_vecs[1][1]/2.0 
        return [xc,yc]
    
    def getUnitCellPos(self,lat): 
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
        return self._get_pos_list(unitCell)
    
    
    
    def _get_pos_list(syst):
        pos_list = []
        for site in syst.sites():
            pos_list.append(site.pos)
        return pos_list 
    
    
    def encode(self,syst):
        pos_list = np.array(self.get_pos_list(syst))
        nlx = round(np.max(pos_list[:,0])/self.lx)
        nly = round(np.max(pos_list[:,1]) - np.min(pos_list[:,1])/self.ly) + 1

        
        
        
        

