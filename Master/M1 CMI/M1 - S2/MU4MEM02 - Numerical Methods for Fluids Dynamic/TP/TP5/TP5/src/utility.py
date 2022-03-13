# coding: utf-8

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg.dsolve import linsolve




def restrict(Afine): 

#=====================================================
#
#  restriction routine (for cell-centered quantities) 
#
#=====================================================

    # resolution
    r  = np.size(Afine,0)
    c  = np.size(Afine,1) 
    n  = r-2;
    m  = c-2; 
    
    n2 = int(n/2);
    m2 = int(m/2);
    
    Acoarse = np.zeros((n2+2,m2+2)) 

    # restriction operation
    for i in range(1,n2+1):
        for j in range(1,m2+1):
            Acoarse[i,j] = (Afine[2*i-2,2*j-2] + \
                            Afine[2*i-1,2*j-2] + \
                            Afine[2*i-2,2*j-1] + \
                            Afine[2*i-1,2*j-1])/4
            
    return Acoarse       



def prolong(Acoarse):

#====================================================
#
#  prolongation routine (for cell-centered quantities) 
#
#====================================================

    # resolution
    n = np.size(Acoarse,0)
    m = np.size(Acoarse,1)
    n2 = 2*(n-2) + 2
    m2 = 2*(m-2) + 2
    s  = (n2,m2)

    Afine = np.zeros(s)

    #prolongation operation
    for i in range(0, n-1):
        for j in range(0, m-1):
            ifine = 2*i-1
            jfine = 2*j-1
            Afine[ifine,jfine] = (9/16*Acoarse[i,j] + \
                                  3/16*Acoarse[i+1,j] + \
                                  3/16*Acoarse[i,j+1] + \
                                  1/16*Acoarse[i+1,j+1])
                            
            Afine[ifine+1,jfine] = (3/16*Acoarse[i,j] + \
                                    9/16*Acoarse[i+1,j] + \
                                    3/16*Acoarse[i+1,j+1] +\
                                    1/16*Acoarse[i,j+1])
        
                              
            Afine[ifine,jfine+1] = (3/16*Acoarse[i,j] + \
                                    1/16*Acoarse[i+1,j] + \
                                    9/16*Acoarse[i,j+1] +
                                    3/16*Acoarse[i+1,j+1]  )
                              
                    
            Afine[ifine+1,jfine+1] = (1/16*Acoarse[i,j] + \
                                      3/16*Acoarse[i+1,j] + \
                                      3/16*Acoarse[i,j+1] + \
                                      9/16*Acoarse[i+1,j+1])
    return Afine
                                 


def diverg(u, v, dt):

    #===========================================
    #
    #  routine computes the divergence
    #  (rhs of pressure Poisson equation) 
    #
    #===========================================

    # resolution
    n = np.size(u,0)
    m = np.size(u,0)
    dx = 1/(n-2)
    dy = 1/(m-2)

    # x-gradient of u
    ux = np.zeros((n,m))
    ux[1:-1,:] = ( u[2:,:] - u[:-2,:] )/(2*dx);

    # y-gradient of v
    vy = np.zeros((n,m))
    vy[:,1:-1] = ( v[:,2:] - v[:,:-2] )/(2*dy);

    # divergence/time-step
    d = np.zeros((n,m))
    d = (ux + vy)/dt
    
    return d


def ucorr(u, p, dt):

    #=============================================
    #
    # routine corrects the u-velocity 
    # (using the x-pressure gradient)
    # 
    #=============================================

    # resolution
    n = np.size(u,0)
    m = np.size(u,1)
    dx = 1/(n-1)

    # compute x-pressure gradient
    px = np.zeros((n,m))
    px[1:-1,:] = ( p[2:,:] - p[:-2,:] )/(2*dx)

    # correction step
    uc = np.zeros((n,m))
    uc = u - dt*px
    
    return uc
    
    
    
def vcorr(v, p, dt):

    #=============================================
    #
    # routine corrects the v-velocity 
    #(using the y-pressure gradient)
    # 
    #=============================================

    # resolution
    n = np.size(v,0)
    m = np.size(v,1)
    dy = 1/(m-2)

    # compute the y-pressure gradient
    py = np.zeros((n,m))
    py[:,1:-1] = ( p[:,2:] - p[:,:-2] )/(2*dy)

    # correction step 
    vc = np.zeros((n,m))
    vc = v - dt*py
    
    return vc




       