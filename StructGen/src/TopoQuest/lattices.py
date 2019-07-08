#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 30 20:33:06 2019

@author: hyda
"""
import kwant


graphene = kwant.lattice.general([[np.sqrt(3)/2,1/2],[np.sqrt(3)/2,-1/2]],  #Lattice vectors 
                                 [[0,0],[1/np.sqrt(3),0]],norbs=1) # Co-ordinates

Zigzag = kwant.lattice.general([[np.sqrt(3)/3,0],[0,1]], #Lattice vectors
                                 [[0,1/6],[np.sqrt(3)/2,2/6],[np.sqrt(3)/2,4/6],[0,5/6]],norbs=1)

Armchair = kwant.lattice.general([[1,0],[0,np.sqrt(3)/3.0]], #Lattice vectors
                                 [[1/6,0],[2/6,np.sqrt(3)/2],[4/6,np.sqrt(3)/2],[5/6,0]],norbs=1)

Lieb = kwant.lattice.general([[1,0],[0,1]],  #Lattice vectors
                             [[0,0],[0.5,0],[0,0.5]],norbs=1) # Coordinates of the sites 

Kagome = kwant.lattice.kagome(norbs=1)
