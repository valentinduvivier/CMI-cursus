# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 15:32:11 2021

@author: sergent
"""

import numpy as np
from scipy.special import erf
from math import pi

def reg_grid(m, xstart, xend): 
    
    """
    s      : regular x distribution
    n      : cell number
    xstart : first point
    xend   : last point
    """
    s  = np.zeros((m))
    
    dh = (xend-xstart)/(m) # constant space step
    s  = np.linspace(xstart, xend, m+1) #np.arange(xstart, xend+dh, dh)
    
    return s

def coordinates(m, x):
    
    """
    m       : cell number
    dx      : scalar cell size
    xp      : scalar node position
    dim_sca : number of scalar nodes
    dxp     : velocity cell size
    xu      : velocity node position
    dim_U   : number of velocity nodes
    dxu     : scalar cell size
    """
    
    # space interval between face positions 
    dx = np.diff(x)
    
    # define scalar nodes inside the domain
    xp = x[:m]+dx/2
    dim_sca = np.array([np.size(xp), 1])
    
    dxp = np.diff(xp)
    
    # ---------------------------------------------
    
    # define velocity nodes inside the domain
    xu = x[1:-1]
    dim_U = np.array([np.size(xu), 1])
    
    dxu = np.diff(xu)
    
    return dx, xp, dim_sca, dxp, xu, dim_U, dxu
