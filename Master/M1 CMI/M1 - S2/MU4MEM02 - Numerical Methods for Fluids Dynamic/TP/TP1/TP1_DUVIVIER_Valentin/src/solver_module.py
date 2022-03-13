#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import scipy.sparse as sp

#-----------------------------------------------------------------

# Gradient and divergence operators
def gradiv(m, dxp):
    
    # declaration
    Div = np.zeros(np.array([m-2,m-1]))
    Grad = np.zeros(np.array([m-1,m]))
    
    # divergence
    Div  = sp.diags([-1, 1], [0, 1], (m-2, m-1)).toarray()
    
    # gradient
    Grad = sp.diags([-1/dxp, 1/dxp], [0, 1], (m-1, m)).toarray()
    
    return Div, Grad

#----------------------------------------------------------------- 

# Interpolation centred scheme - face centred mesh
def interp(m):

   # declaration
   Int = np.zeros(np.array([m-1, m]))

   # interpolation
   Int = sp.diags([1/2, 1/2], [0, 1], (m-1, m)).toarray()
   
   return Int







