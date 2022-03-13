#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# BOUNDARY CONDITIONS
"""

import numpy as np
import scipy.sparse as sp

#-------------------------------------------------------------------------------

# diffusion for Phi at VF interface
def bc_diff(bcW, bcE, SbcW, SbcE, m, x, xp, dxp, C0, CL, gamma_x, Gamx_ew, mass):
       
    D  = gamma_x/dxp[0]
    
    aW = 0
    aE = D
    
    SP = -(2*D)
    Su = 2*D * C0
    
    # -------------------
    
    ## Combination
    aP = (aW + aE - SP)
    bP = Su
    
    # -------------------
    
    # WEST node
    bcW[0,0] = - aP
    SbcW = - bP
    
    # -------------------
    # -------------------
    
    aW = D
    aE = 0
    
    SP = -(2*D)
    Su = 2*D * CL
    
    # -------------------
    
    ## Combination
    aP = (aW + aE - SP)
    bP = Su
    
    # -------------------
    
    # EST node
    bcE[0,-1] = - aP
    SbcE = - bP
    
    return bcW, bcE, SbcW, SbcE

#-------------------------------------------------------------------------------

# centred interpolation for Phi at VF interface
def bc_interp(bcW, bcE, SbcW, SbcE, U0, UL, C0, CL, rhoval, rho_ew, Ux, mass):
       
    F  = rhoval*Ux[0]
    
    aW = 0
    aE = -F/2
    
    SP = -F
    Su = F*C0
    
    # -------------------
    
    ## Combination
    aP = (aW + aE - SP)
    bP = Su
    
    # -------------------
    
    # WEST node
    bcW[0,0] += - aP
    SbcW += - bP
    
    # -------------------
    
    F = rhoval*Ux[0]

    aW = F/2
    aE = 0
    
    SP = -(-F)
    Su = (-F)*CL
    
    # -------------------

    ## Combination
    aP = (aW + aE - SP)
    bP = Su

    # -------------------    
    # EST node
    bcE[0,-1] += - aP
    SbcE += - bP
    
    return bcW, bcE, SbcW, SbcE

#-------------------------------------------------------------------------------



    