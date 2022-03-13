# coding: utf-8

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg.dsolve import linsolve
from utility import *
 

def Poiss(p_old, f):
    
    # resolution
    n = np.size(f,0)
    m = np.size(f,1)
    
    s = (n,m)
    
    p       = np.zeros(s)
    p_new   = np.zeros(s)
    
    # number of multigrid cycles
    ncycle = 10
    
    p = p_old   # initialize with previous pressure field 
    
    for i in range(1,ncycle+1): 
        p = MGV_Poi(p, f)
        p = bc_neu(p)
        r = residual_Poi(p, f)
        
        # scale the righthand side (necessary because of Neumann data on the edge) 
        # this scaling makes the solution unique
        f[1:-1,1:-1] = f[1:-1,1:-1] - np.sum(r[1:-1,1:-1])/(n-2)/(m-2);
    
    p_new = bc_neu(p)   
    
    return p_new
    

    
def MGV_Poi(Ain, f):

    #=====================================================
    #
    # geometric multigrid V-cycle for Poisson equation
    #
    # (recursive definition) 
    # (dimension N has to be of the form 2^k + 2) 
    #
    # Ain:  guessed solution (n x m)-matrix
    # f:    right-hand side (n x m)-matrix
    #
    #=====================================================
    
    # resolution
    n = np.size(f,0); 
    m = np.size(f,1);


    # if we are at the coarsest level
    if ((n==4)and(m==4)) :
        Aout = relax_Poi(Ain,f,10) 
    else : 
    # otherwise

        # relax 10 times (pre-smoothing)
        Asmooth = relax_Poi(Ain, f, 10)

        #compute the residual
        res  = residual_Poi(Asmooth, f) 

        #restrict the residual to the next-coarser grid
        res2 = restrict(res) 

        #solve the error equation on the next-coarser grid by MGV
        err = MGV_Poi(np.zeros(np.shape(res2)), res2) 
        
        #add the prolongated error to the solution 
        Asmooth = Asmooth + prolong(err)

        #relax 10 times (post-smoothing) 
        Aout = relax_Poi(Asmooth, f, 10) 

    return Aout


    
def relax_Poi(xguess, b, ktimes):

#===========================================
#
# solves the pressure Poisson equation 
# using Gauss-Seidel
#
#===========================================

    # resolution
    n     = np.size(b,0)
    m     = np.size(b,1)
    dx = 1/(n-2)
    dy = 1/(m-2)

    coefx = 1/dx/dx
    coefy = 1/dy/dy
    coef0 = 2*(coefx + coefy)
    
    # initialization
    p = np.zeros((n,m))
    p = xguess;

    # iteration
    for k in range(1, ktimes+1):
        for i in range(1, n-1): 
            for j in range(1, m-1):
                p[i,j] = (coefx*(p[i+1,j]+p[i-1,j]) + \
                          coefy*(p[i,j+1]+p[i,j-1]) - \
                          b[i,j])/coef0
        p = bc_neu(p)
      
    return p
    
    
def residual_Poi(Ain, f) :   
#================================================
#
#  residual routine for pressure Poisson equation 
#
#================================================
    # resolution
    n     = np.size(f,0)
    m     = np.size(f,1)
    
    dx    = 1/(n-2)
    dy    = 1/(m-2)

    # coefficients for the Poisson equation
    coefx = 1/dx/dx
    coefy = 1/dy/dy
    coef0 = 2*(coefx + coefy)

    s     = (n,m)
    res   = np.zeros(s)
    

    # implement boundary conditions
    Ain = bc_neu(Ain)
    
    
    # residual computation
    res = np.zeros(np.shape(Ain)); 
    res[1:-1,1:-1] = (f[1:-1,1:-1] + Ain[1:-1,1:-1]*coef0 - \
                      (Ain[1:-1,2:] + Ain[1:-1,:-2])*coefx - \
                      (Ain[2:,1:-1] + Ain[:-2,1:-1])*coefy )
    
    
    return res


def bc_neu(q):

    #------------------------------------------
    # implement homogeneous Neumann boundary conditions
    #------------------------------------------

    # ghost cell mapping
    q[:, 0]  = q[:,  1]
    q[:,-1]  = q[:, -2]
    q[0, :]  = q[1,  :]
    q[-1,:]  = q[-2, :]

    # corner elements (only needed for prolongation)
    q[0,  0] = q[1,  1] 
    q[-1, 0] = q[-2, 1]
    q[0, -1] = q[1, -2]
    q[-1,-1] = q[-2,-2]
     
        
    return q    
     