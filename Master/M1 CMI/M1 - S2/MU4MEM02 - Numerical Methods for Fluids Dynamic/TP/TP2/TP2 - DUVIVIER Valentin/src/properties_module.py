#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

def prop_phys(dim_sca, dim_U, gamma_x, rhoval):
    
    """
    Gamx_ew : diffusion coeffient on the faces
    rho_ew  : density on the faces
    """
    
    # diffusion coeffient at the scalar nodes
    Gamx = gamma_x * np.ones(dim_sca)
    rho  = rhoval  * np.ones(dim_sca)
    
    # diffusion coeffient on the faces
    Gamx_ew = np.zeros(dim_U) # harmonic average
    rho_ew  = np.zeros(dim_U) # arithmetic average
    
    if gamma_x > 0:
        Gamx_ew = 2 * (Gamx [0:-1,:]*Gamx [1:,:]) / (Gamx [0:-1,:] + Gamx[1:,:])

    if rhoval > 0:
        rho_ew = .5 * (rho [0:-1,:] + rho [1:,:]) 
    
    return Gamx_ew, rho_ew
