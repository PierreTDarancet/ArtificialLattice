#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 16:27:19 2019

@author: ssrinivasan
"""

#from helper import *
from helper import Armchair,get_width,finite_to_1D
import kwant
import random 
import numpy as np 

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
    
    def _get_random_2pts(self,L,w): 
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
    
    def make_symmetric(self,symmetry=['mirror'],Ncentral=7): 
               
        min_width = get_width(Ncentral,self.lat)
        rect_L_plus_W = min_width + self.lx
        offset =self.ly
        
        def _shape_from_lines(pos,lineCoeff1,lineCoeff2,offset=offset):
                x,y = pos
                val1 = np.polyval(lineCoeff1,abs(x-self.lx))
                val2 = np.polyval(lineCoeff2,abs(x-self.lx))
                if 0.0<=y+offset<self.ly and  0<=x<2*self.lx:
                    if val1<=y<=val2: 
                        return True 
                    else: 
                        return False
                else: 
                        return False
                
        while offset > self.ly-self.lat.prim_vecs[1][1]: 
            ypt11,ypt12 = self._get_random_2pts(rect_L_plus_W,min_width)
            ypt11,ypt12 = min(ypt11,ypt12),max(ypt11,ypt12)
            PTS1 = []
            for pt in [ypt11,ypt12]: 
                if pt > self.ly:
                    PTS1.append([pt-self.ly,self.ly])
                else: 
                    PTS1.append([0,pt])
                
                ypt21,ypt22 = self._get_random_2pts(self.ly,min_width)
                ypt21,ypt22 = min(ypt21,ypt22),max(ypt21,ypt22)
                PTS2= [[self.lx,ypt21],[self.lx,ypt22]]
        
        
            lineCoeff1 =  np.polyfit([PTS1[0][0],PTS2[0][0]],
                                     [PTS1[0][1],PTS2[0][1]],1)
        
            lineCoeff2 =  np.polyfit([PTS1[1][0],PTS2[1][0]],
                                     [PTS1[1][1],PTS2[1][1]],1)
            offset = (lineCoeff1[1]+lineCoeff2[1])/2
            lineCoeff1[1] -= offset 
            lineCoeff2[1] -= offset
        
        syst = kwant.Builder(kwant.TranslationalSymmetry([2*self.lx,0]))
        syst[self.lat.shape((lambda pos: _shape_from_lines(pos,
                                              lineCoeff1=lineCoeff1, 
                                              lineCoeff2=lineCoeff2,
                                              offset=offset)),(0.0,0))]=0 
        syst[self.lat.neighbors()]=-1
        #syst = finite_to_1D(syst,self.lx)
        #syst.eradicate_dangling()
        return syst
            
            
        

        