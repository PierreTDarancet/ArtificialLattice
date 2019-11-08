#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 17:04:37 2019

@author: ssrinivasan
"""

""" Module to train logistic regression to predict the topological invariant
    from Motif Matrix 
"""
from sklearn.linear_model import LogisticRegression
import numpy as np 

class logistic(object):
    def __init__(self,gridCV=False,solver='liblinear',C=1.0,penalty='l2'):
        self.reg=LogisticRegression(solver=solver,
                                    C=C,
                                    penalty=penalty)
        if gridCV is False:
            self.fit = self.reg.fit
        else: 
            from sklearn.model_selection import GridSearchCv 
            parameters = {'C':np.lo(0.1,100,5)}
            clf = GridSearchCv(self.reg,parameters,cv=5)
            self.fit=clf.fit
            

        