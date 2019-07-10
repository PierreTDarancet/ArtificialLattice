#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 18:24:52 2019

@author: ssrinivasan
"""

class sectioned(general):
    """
    For reading file type with sections. The subroutine 
    records byte position of each line and the line 
    number that mark the beginning of a section.
    """

    def __init__(self,fname,**kwargs):
        """
        only line < Nmax are potential headline 
        (for efficiency purpose)

        self.
        line      start byte position of each line
        sect      start and end line number of all sections
        sections  all available sections
        """
        # process args
        self.mode = kwargs.get('mode','r')
        self.f = open(fname,self.mode)
        self.fname = fname
        # initialize
        self.header = kwargs.get('headder','ITEM: TIMESTEP')
        self.bytes = [[0,None]]
#        self.line, self.sect = [0],{'header':[1,None]}
        self._parse()


    def _parse(self):
        """
        Parse byte positions of linebreak.
        This is used for seeking file with variable 
        line length.

        This procedure also process sections
        """
        Nbyte = 0 # byte counter
        for line in self.f:
            # byte position of linebreak
            Nbyte += len(line)
            self.line.append(Nbyte)
            self.bytes.append([Nbyte,None])
            # last line of last section = end of file
            self.bytes[-2][1] = Nbyte
        # rewind the file
        self.f.seek(0)
        return

    def read(self,frame):
        """
        Extract section
        use .decode("utf-8") to get nice printing
        """
        start = self.bytes[frame][0]
        end   = self.line[self.sect[section][1]]-1
        self.f.seek(start)