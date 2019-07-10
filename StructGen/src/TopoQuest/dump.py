#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 16:04:52 2019

@author: ssrinivasan
"""

class dump(): 
    def __init__(self,filename='dump.atom'): 
        self.filename = filename 
        self.f = open(self.filename,'w+')
        
    def write_frame(self,timestep,lammpsBox,natoms,atoms_data):
        
        dump_string = ['ITEM: TIMESTEP'] 
        dump_string.append('{}'.format(timestep))
        dump_string.append('ITEM: NUMBER OF ATOMS')
        dump_string.append('{}'.format(natoms))
        dump_string.append('ITEM: BOX BOUNDS pp pp pp')
        box_string = ['{} {}'.format(lammpsBox[0][0],lammpsBox[0][1]), 
                      '{} {}'.format(lammpsBox[1][0],lammpsBox[1][1]),
                      '{} {}'.format(lammpsBox[2][0],lammpsBox[2][1])]
        dump_string += box_string
        attributes = atoms_data['attributes']
        dump_string.append('ITEM: ATOMS ' + " ".join(attributes))
        atom_string = []
        adata = atoms_data['data']
        for i in range(len(adata)):
            atom_string.append(' '.join(str(item) for item in adata[i]))
        dump_string += atom_string 
        self.f.write("\n".join(dump_string))
        self.f.write("\n")
        

            