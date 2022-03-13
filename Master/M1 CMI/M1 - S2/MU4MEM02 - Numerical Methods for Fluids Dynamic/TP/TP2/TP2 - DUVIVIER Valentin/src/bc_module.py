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
    
    # WEST node
    bcW[0,0] = bcW[0,0] + mass[0,0] * gamma_x / (xp[0]-x[0]) - mass[0,0] *Gamx_ew[0] / dxp[0]
    SbcW = SbcW + mass[0,0] * gamma_x / (xp[0]-x[0]) * C0
    
    # EST node
    bcE[0,-1] = bcE[0,-1] + mass[0,0] *gamma_x / (x[-1]-xp[-1]) - mass[0,0] *Gamx_ew[-1] / dxp[-1] 
    SbcE = SbcE - mass[0,0] * gamma_x / (x[-1]-xp[-1]) * CL
    
    return bcW, bcE, SbcW, SbcE

#-------------------------------------------------------------------------------

# diffusion quick for Phi at VF interface
def bc_diff_quick(bcW0, bcEN_1,  SbcW0, SbcEN_1, m, x, xp, dxp, C0, CL, gamma_x, Gamx_ew, mass, dx):
    
    # WEST node
    
    # condition in i=0
    bcW0[0,0] = bcW0[0,0] - (- gamma_x/dxp[0] * xp[1]/xp[0] - gamma_x/dxp[0])
    bcW0[0,1] = bcW0[0,1] - (  gamma_x/dxp[0] * xp[0]/xp[1] + gamma_x/dxp[0])

    SbcW0     = SbcW0 - gamma_x*(dxp[0] + dxp[1])/(xp[0]*xp[1]) * C0

    # ---------------------------------------------------------------------------     
    # EST node
    
    # condition in i=n-1
    bcEN_1[0,-1] = bcEN_1[0,-1] + (gamma_x*((dxp[-1] + dx[-1]/2) / (dxp[-1] * dx[-1]/2)) + gamma_x/dxp[-1])
    bcEN_1[0,-2] = bcEN_1[0,-2] - (gamma_x*(dx[-2]/2 / (dxp[-1] * (dxp[-1]+dx[-1]/2))) + gamma_x/dxp[-1])

    SbcEN_1     = SbcEN_1 + (gamma_x*((dxp[-1] + dx[-1]) / (dx[-1]/2 * (dxp[-1]+dx[-1]/2)))) * CL

    return bcW0, bcEN_1,  SbcW0, SbcEN_1

# -------------------------------------------------------------------------------

# centred interpolation for Phi at VF interface
def bc_interp(bcW, bcE,SbcW, SbcE, U0, UL, C0, CL, rhoval, rho_ew, Ux, mass):
    
    # WEST node
    bcW[0,0] = bcW[0,0] + mass[0,0] * (rho_ew*Ux)[0]/2.
    SbcW = SbcW + mass[0,0]   * rhoval*U0 * C0 
    
    # EST node
    bcE[0,-1] = bcE[0,-1] - mass[-1,-1] * (rho_ew*Ux)[-1]/2. 
    SbcE = SbcE - mass[-1,-1] * rhoval*UL * CL 
    
    return bcW, bcE, SbcW, SbcE

#-------------------------------------------------------------------------------

# upwind scheme for for Phi at VF interface
def bc_upwind(bcW, bcE, SbcW, SbcE, U0, UL, C0, CL, rhoval, rho_ew, Ux, mass):
    
    # WEST node
    
    bcW[0,0] = bcW[0,0] - (mass[1,1] * (rho_ew*Ux)[1] * (Ux[1]>0) - mass[0,0]*(rho_ew*Ux)[0] * (Ux[0]<0))  +  mass[1,1] * (Ux[1]>0)*(rho_ew*Ux)[1]
    
    SbcW = SbcW + mass[0,0]*rhoval * U0 * C0
    
    # ---------------------------------------------------------------------------    
    # EST node

    bcE[0,-1] = bcE[0,-1] - (mass[-2,-2] * (rho_ew*Ux)[-2] * (Ux[-2]>0) - mass[-1,-1]*(rho_ew*Ux)[-1] * (Ux[-1]<0))  +  mass[-2,-2] * (Ux[-2]>0)*(rho_ew*Ux)[-2]
    
    SbcE = SbcE - mass[-1,-1]*rhoval * UL * CL
    
    return bcW, bcE, SbcW, SbcE

#-------------------------------------------------------------------------------

def bc_hybrid(bcW, bcE, SbcW, SbcE, U0, UL, C0, CL, rhoval, rho_ew, Ux, mass, dx, Gamx_ew):
    
    # WEST node

    RW = rho_ew[0]*U0*dx[0] / Gamx_ew[0]
    
    if np.abs(RW) < 2:
        bcW[0,0] = bcW[0,0] + mass[0,0] * (rho_ew*Ux)[0]/2.
        SbcW     = SbcW     + mass[0,0] * rhoval*U0 * C0 
    
    else:
        bcW[0,0] =  mass[0,0] * rho_ew[0] * Ux[0] * (Ux[0]>0) + 2*(Gamx_ew[0]/dx[0])

        SbcW     = (mass[0,0] * rhoval*U0  + 2*(Gamx_ew[0]/dx[0])) * C0

    # ---------------------------------------------------------------------------     
    # EST node
    
    RE = rho_ew[-1]*U0*dx[-1] / Gamx_ew[-1]

    if np.abs(RE) < 2:
        bcE[0,-1] = bcE[0,-1] - mass[-1,-1] * (rho_ew*Ux)[-1]/2. 
        SbcE      = SbcE      - mass[-1,-1] * rhoval*UL * CL 
    
    else:
        bcE[0,-1] = mass[-1,-1] * rho_ew[-1] * Ux[-1] * (Ux[-1]>0) + 2*(Gamx_ew[-1]/dx[-1])

        SbcE      = (mass[-1,-1] * rhoval*UL + 2*(Gamx_ew[-1]/dx[-1])) * CL
        
    return bcW, bcE, SbcW, SbcE

#-------------------------------------------------------------------------------

# quick scheme for for Phi at VF interface
def bc_quick(bcW0, bcW1, bcEN_1, bcEN_2,  SbcW0, SbcW1, SbcEN_1, SbcEN_2, U0, UL, C0, CL, rhoval, rho_ew, Ux, mass):
    
    # WEST node
    
    # condition in i=0
    bcW0[0,0] = bcW0[0,0] - (-((Ux[0] > 0)* 7 + (Ux[0] < 0)* 3)) * rho_ew[0] * Ux[0] / 8
    bcW0[0,1] = bcW0[0,1] + (  (Ux[0] > 0)* 3 + (Ux[0] < 0)* 6)  * rho_ew[0] * Ux[0] / 8

    SbcW0     = SbcW0 + (- mass[0,0] * rhoval*U0 - mass[0,0] * rho_ew[0] * Ux[0] * (Ux[0] > 0) / 4)* C0 

    # condition in i=1
    bcW1[0,0] = bcW1[0,0] + (- 7 * rho_ew[0] * Ux[0] * (Ux[0] > 0) / 8 - rho_ew[1] * Ux[1] * (Ux[1] > 0) / 8)
    
    SbcW1     = SbcW1 + mass[1,1] * rho_ew[0]*Ux[0]* C0 / 4 * (Ux[0] > 0)
    
    # ---------------------------------------------------------------------------     
    # EST node
    
    # condition in i=n-1
    bcEN_1[0,-1] = bcEN_1[0,-1] + (-((Ux[-1] > 0)* 3 + (Ux[-1] < 0)* 7)) * rho_ew[-1] * Ux[-1] / 8
    bcEN_1[0,-2] = bcEN_1[0,-2] - (  (Ux[-2] < 0)* 3 + (Ux[-2] > 0)* 6)  * rho_ew[-2] * Ux[-2] / 8
    bcEN_1[0,-3] = bcEN_1[0,-3] + (  (Ux[-1] > 0)  * rho_ew[-1] * Ux[-1] / 8)
    
    SbcEN_1      = SbcEN_1 + (mass[-1,-1] * rhoval*UL + mass[-2,-2] * rho_ew[-2] * Ux[-2] * (Ux[0] < 0) / 4)* CL

    # condition in i=n-2
    bcEN_2[0,-1] = bcEN_2[0,-1] + (Ux[0] < 0)* 7 * mass[-2,-2] * rho_ew[-1] * Ux[-1] / 8 + mass[-2,-2] * rho_ew[-2] * Ux[-2] * (Ux[-2] < 0) / 8

    SbcEN_2      = SbcEN_2 + (-mass[1,1] * rho_ew[-1]*Ux[-1]* CL * (Ux[-2] < 0)/ 4)
        
    return bcW0, bcW1, bcEN_1, bcEN_2,  SbcW0, SbcW1, SbcEN_1, SbcEN_2

