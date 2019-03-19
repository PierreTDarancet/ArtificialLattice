#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 16:27:19 2019

@author: ssrinivasan
"""

from helper import *
import kwant
from random

class GenStruct(): 
    def __init__(self): 
        self.lat = Armchair
        self.lx = 20*Armchair.prim_vecs[0][0]
        self.ly = 20*Armchair.prim_vecs[1][1]
        self.onsite = 0
        self.hop = -1
        self.nsites = 500
        self.syst = None
            
    def make_finite(self): 
        syst = kwant.Builder()
        self.make_symmetric(syst)
        return syst
    
    def make_symmetric(self,syst,symmetry=['inversion'],Ncentral=5):
        
        Ncentral = 5
        ymin = get_width(Ncentral,self.lat)
        x1,y1 = random.uniform(0,lx),random.uniform(ymin,ly)
        x2,y2 = random.unifrm
        
        
        def fill_central_strip(site):
            x,y = site.pos
            ymin = self.lat.prim_vecs[1][1]
            if 0<= x < self.lx/2 and 0<=y< ymin: 
                return True
            else: 
                return False 
        syst[self.lat(fill_central_strip,(0,0))] = self.onsite 
        n_filled = len(list(syst.sites()))
        n_maxfill = self.nsites/2.0
        while nfilled < n_maxfill: 
            site_nucleate = random.choice(list(syst.sites()))
            
            
        

        