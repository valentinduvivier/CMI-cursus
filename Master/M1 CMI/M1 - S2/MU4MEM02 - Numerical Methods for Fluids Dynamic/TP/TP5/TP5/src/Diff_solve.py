# coding: utf-8

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg.dsolve import linsolve
from utility import *
 

def DiffU(uold, p, dt, n):

    #=======================================
    # solve a diffusion problem for the 
    # new velocity (implicit; backward Euler)
    #=======================================
    
    # resolution
    n = np.size(uold,0)
    m = np.size(uold,1)
    
    s = (n,m)
    
    u      = np.zeros(s)
    unew   = np.zeros(s)

    ncycle = 6                  # number of multigrid cycles 

    u      = uold               # use previous field as initial guess 

    for i in range(1,ncycle+1):
        u = MGV_Diff_U(u, uold, p, dt, n)

    unew = bc_diri_U(u, p, dt, 1)  # impose boundary conditions 
                     
    return unew	


def DiffV(vold, p, dt, n): 

    #=======================================
    #
    # solve a diffusion problem for the 
    # new velocity (implicit; backward Euler)
    #
    #=======================================
    
    # resolution
    n = np.size(vold, 0)
    m = np.size(vold, 1)
    
    s = (n, m)
    
    v      = np.zeros(s)
    vnew   = np.zeros(s)

    ncycle = 6                      # number of multigrid cycles

    v      = vold                   # use previous field as initial guess 

    for i in range(1, ncycle+1):
        v = MGV_Diff_V(v, vold, p, dt, n)


    vnew = bc_diri_V(v, p, dt, 1)     # impose boundary conditions
 
    return vnew



def MGV_Diff_U(Ain, f, p, dt, nmax) :

    #------------------------------------------------------------------
    # geometric multigrid V-cycle for Diffusion equation
    #
    # (recursive definition) 
    # (dimension N has to be of the form 2^k + 2) 
    #
    # Ain:  guessed solution (n x m)-matrix
    # f:    right-hand side (n x m)-matrix
    # p:    pressure (needed for boundary conditions) 
    # dt:   time-step 
    # nmax: resolution at maximum level
    #
    #------------------------------------------------------------------

    n = np.size(f, 0) 
    m = np.size(f, 1)  
     
    if (n == nmax):      # set the boundary condition flag 
        bcflag = 1       # for top level equations (inhomogeneous)
        
    else:
        bcflag = 0       # for error equations (homogeneous)

    # if we are at the coarsest level
    if ((n==4)and(m==4)):  
        Aout = relax_Diff_U(Ain, f, p, dt, 10, bcflag)
        
    else:
    # otherwise

        #relax 10 times (pre-smoothing)
        Asmooth = relax_Diff_U(Ain, f, p, dt, 10, bcflag) 

        #compute the residual
        res     = residual_Diff_U(Asmooth, f, p, dt, bcflag)

        #restrict the residual and the pressure to the next-coarser grid
        res2    = restrict(res) 
        p2      = restrict(p)
        res2    = bc_diri_U(res2, p2, dt, 0)

        #solve the error equation on the next-coarser grid by MGV
        err     = MGV_Diff_U(np.zeros(np.shape(res2)), res2, p2, dt, nmax)

        #add the prolongated error to the solution 
        err     = bc_diri_U(err, p2, dt, 0)
        Asmooth = Asmooth + prolong(err) 

        #relax 10 times (post-smoothing) 
        Aout    = relax_Diff_U(Asmooth, f, p, dt, 10, bcflag)
        
    #impose boundary conditions
    Aout = bc_diri_U(Aout, p, dt, bcflag)
    
    return Aout



def MGV_Diff_V(Ain, f, p, dt, nmax):

    #------------------------------------------------------------------
    # geometric multigrid V-cycle for Diffusion equation
    #
    # (recursive definition) 
    # (dimension N has to be of the form 2^k + 2) 
    #
    # Ain:  guessed solution (n x m)-matrix
    # f:    right-hand side (n x m)-matrix
    # p:    pressure (needed for boundary conditions) 
    # dt:   time-step 
    # nmax: resolution at maximum level
    #
    #------------------------------------------------------------------

    n = np.size(f,0) 
    m = np.size(f,1)  

    if (n == nmax):   # set the boundary condition flag 
        bcflag = 1    # for top level equations (inhomogeneous)
    else:
        bcflag = 0    # for error equations (homogeneous) 

    # if we are at the coarsest level
    if ((n==4)and(m==4)): 
        Aout = relax_Diff_V(Ain, f, p, dt, 10, bcflag) 
    
    else:
    # otherwise

        #relax 10 times (pre-smoothing)
        Asmooth = relax_Diff_V(Ain, f, p, dt, 10, bcflag) 

        #compute the residual
        res  = residual_Diff_V(Asmooth, f, p, dt, bcflag) 

        #restrict the residual and the pressure to the next-coarser grid
        res2 = restrict(res) 
        p2   = restrict(p)
        res2 = bc_diri_V(res2,p2,dt,0)

        #solve the error equation on the next-coarser grid by MGV
        err = MGV_Diff_V(np.zeros(np.shape(res2)), res2, p2, dt, nmax) 

        #add the prolongated error to the solution 
        err = bc_diri_V(err,p2,dt,0)
        Asmooth = Asmooth + prolong(err) 

        #relax 10 times (post-smoothing) 
        Aout = relax_Diff_V(Asmooth, f, p, dt, 10, bcflag)


    # impose boundary conditions
    Aout = bc_diri_V(Aout, p, dt, bcflag)
    
    return Aout



def relax_Diff_U(xguess, b, p, dt, ktimes, bcflag):

    #=============================================
    #
    #    solves the diffusion equation for u-velocity 
    #    using Gauss-Seidel
    #
    #=============================================

    # resolution
    n = np.size(b,0)
    m = np.size(b,1)
    
    dx = 1/(n-2)
    dy = 1/(m-2)

    # coefficients for diffusion equation
    coefx = dt/dx/dx
    coefy = dt/dy/dy
    coef0 = 1 + 2*(coefx + coefy) 

    # initialization
    u = np.zeros(np.shape(xguess))
    u = xguess

    # implement boundary conditons
    u = bc_diri_U(u, p, dt, bcflag)

    # iteration
    
    for k in range(1,ktimes+1):
        for i in range(1,n-1):
            for j in range (1,m-1):
                u[i,j] = (coefx*(u[i+1,j]+u[i-1,j]) + \
                          coefy*(u[i,j+1]+u[i,j-1]) + \
                          b[i,j])/coef0
     
        u = bc_diri_U(u, p, dt, bcflag)
   

    return u


def relax_Diff_V(xguess, b, p, dt, ktimes, bcflag):

    #=============================================
    #
    #    solves the diffusion equation for v-velocity 
    #    using Gauss-Seidel
    #
    #=============================================

    # resolution
    n = np.size(b,0)
    m = np.size(b,1)
    dx = 1/(n-2)
    dy = 1/(m-2)

    # coefficients for diffusion equation
    coefx = dt/dx/dx
    coefy = dt/dy/dy
    coef0 = 1 + 2*(coefx + coefy) 

    # initialization
    v = np.zeros(np.shape(xguess))
    v = xguess

    # implement boundary conditons
    v = bc_diri_V(v,p,dt,bcflag)

    for k in range(1,ktimes+1):
        for i in range(1,n-1):
            for j in range (1,m-1):
                v[i,j] = (coefx*(v[i+1,j]+v[i-1,j]) + \
                          coefy*(v[i,j+1]+v[i,j-1]) + \
                          b[i,j])/coef0
     
        v = bc_diri_V(v, p, dt, bcflag)
        
    return v


def residual_Diff_U(Ain, f, p, dt, bcflag)  :

    #================================================
    #
    #  residual routine for diffusion equation 
    #  (u-velocity)
    #
    #================================================

    # resolution
    n = np.size(f,0)
    m = np.size(f,1)
    dx = 1/(n-2)
    dy = 1/(m-2)

    # coefficients for the diffusion equation
    coefx = dt/dx/dx
    coefy = dt/dy/dy
    coef0 = 1 + 2*(coefx + coefy) 

    # implement boundary conditions
    Ain = bc_diri_U(Ain, p, dt, bcflag) 

    # residual computation
    res = np.zeros(np.shape(Ain)) 
    res[1:-1,1:-1] = (f[1:-1,1:-1] - Ain[1:-1,1:-1]*coef0 + \
                      (Ain[1:-1,2:] + Ain[1:-1,:-2])*coefx + \
                      (Ain[2:,1:-1] + Ain[:-2,1:-1])*coefy)
    
    return res

    
def residual_Diff_V(Ain, f, p, dt, bcflag)  :

    #================================================
    #
    #  residual routine for diffusion equation 
    #  (v-velocity)
    #
    #================================================

   # resolution
    n = np.size(f, 0)
    m = np.size(f, 1)
    dx = 1/(n-2)
    dy = 1/(m-2)

    # coefficients for the diffusion equation
    coefx = dt/dx/dx
    coefy = dt/dy/dy
    coef0 = 1 + 2*(coefx + coefy) 

    # implement boundary conditions
    Ain = bc_diri_V(Ain, p, dt, bcflag) 

    # residual computation
    res = np.zeros(np.shape(Ain)) 
    res[1:-1,1:-1] = (f[1:-1,1:-1] - Ain[1:-1,1:-1]*coef0 + \
                      (Ain[1:-1,2:] + Ain[1:-1,:-2])*coefx + \
                      (Ain[2:,1:-1] + Ain[:-2,1:-1])*coefy)
    
    return res

def  bc_diri_U(u, p, dt, bcflag): 

    #===========================================
    #
    # implement Dirichlet boundary conditions
    # (for u-velocity with pressure correction)
    #
    #===========================================

    # resolution
    n  = np.size(u,0);
    m  = np.size(u,1);
    dx = 1/(n-2);
    dy = 1/(m-2);

    # compute x-pressure gradient along sides 
    px1 = np.zeros(n)
    px1[1:-1] = (p[2:,0] - p[:-2,0])/(2*dx)
    px2 = np.zeros(n)
    px2[1:-1] = (p[2:,-1] - p[:-2,-1])/(2*dx)     
        
    # ghost cell mapping
    u[0, :] = -u[1, :]   
    u[-1,:] = -u[-2,:]  
    u[:, 0] = -u[:, 1] + 2*bcflag*(      dt*px1)
    u[:,-1] = -u[:,-2] + 2*bcflag*(1.0 + dt*px2)

    # corner elements (only needed for prolongation)
    u[0,  0] = u[1,  1] 
    u[-1, 0] = u[-2, 1]
    u[0, -1] = u[1, -2]
    u[-1,-1] = u[-2,-2]
        
    return u

def bc_diri_V(v, p, dt, bcflag): 

    #===========================================
    #
    # implement Dirichlet boundary conditions
    # (for v-velocity with pressure correction)
    #
    #===========================================

    # resolution
    n  = np.size(v,0)
    m  = np.size(v,0)
    dx = 1/(n-2)
    dy = 1/(m-2)

    # compute y-pressure gradient along side 
    py1 = np.zeros(m)
    py1[1:-1] = (p[0,2:] - p[0,:-2])/(2*dy)
    py2 = np.zeros(m)
    py2[1:-1] = (p[-1,2:] - p[-1,:-2])/(2*dy) 

    # ghost cell mapping
    v[0,  :] = -v[1,  :]  + 2*bcflag*(      dt*py1) 
    v[-1, :] = -v[-2, :]  + 2*bcflag*(1.0 + dt*py2)  
    v[:,  0] = -v[:,  1] 
    v[:, -1] = -v[:, -2] 

    # corner elements (only needed for prolongation)
    v[0,  0] =  v[1,  1] 
    v[-1, 0] =  v[-2, 1]
    v[0, -1] =  v[1, -2]
    v[-1,-1] =  v[-2,-2]
   
    return v
    