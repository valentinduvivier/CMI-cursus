# coding: utf-8

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg.dsolve import linsolve

#===============================================
#
#  routine solves pressure Poisson equation 
#  (with Dirichlet boundary conditions) 
#
#===============================================


def Poiss(phi_exact, phi_old, q, ncycle):
    
    # resolution
    n = np.size(q,0)
    m = np.size(q,1)
    
    s = (n,m)
    
    phi_new = np.zeros(s)
    p       = np.zeros(s)
        
    p = phi_old
    
    # ---------------------------
    
    error    = np.sum(np.abs(phi_old - phi_exact))/np.sum(phi_exact)
    counter  = 0
    
    # ---------------------------
    
    while error > 10**-6:
        p = MGV_Poi(p, q, ncycle)
        r = residual_Poi(p, q)
        
        print("The residual =", r.max())
        
        counter += 1
        error    = np.sum(np.abs(p - phi_exact))/np.sum(phi_exact)
        
        if counter == ncycle:
            print(f'Bad : max iterations reached, counter = {counter}')
            break
    
    print(f'\nGood : error reached, counter = {counter}\n')
    
    phi_new = p
    
    return phi_new
 
    
    
def MGV_Poi(Ain, f, mu):

#=====================================================
#
# geometric multigrid V-cycle for Poisson equation
#
# (recursive definition) 
#
# Ain:  guessed solution (n x m)-matrix
# f:    right-hand side (n x m)-matrix
#
#=====================================================

# resolution
    n = np.size(f, 0); 
    m = np.size(f, 1);
        
    # if we are at the coarsest level
    if ((n==4)and(m==4)) :
        Aout = relax_Poi(Ain, f, mu) 

    else :
    # otherwise
        #relax 10 times (pre-smoothing)
        v = relax_Poi(Ain, f, mu)

        #compute the residual
        r = residual_Poi(v, f)

        #restrict the residual to the next-coarser grid
        R = restrict(r)
        
        #solve the error equation on the next-coarser grid by MGV
        E = MGV_Poi(np.zeros(np.shape(R)), R, mu)
        
        #add the prolongated error to the solution 
        v += prolong(E)

        #relax 10 times (post-smoothing) 
        Aout = relax_Poi(v, f, mu)

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
    delta = 2/(n-1)

    # initialization
    p = xguess;

    # iteration
    for k in range(1, ktimes+1):
        for i in range(1,n-1): 
            for j in range(1,m-1):
                p[i,j] = (1/4)*(p[i-1,j] + p[i+1,j] + p[i,j-1] + p[i,j+1]) + (delta**2/4)*(b[i,j])

    # final solution
    x = p
      
    return x
    
    
def residual_Poi(Ain, f) :   
#================================================
#
#  residual routine for pressure Poisson equation 
#
#================================================
    # resolution
    n     = np.size(f,0)
    m     = np.size(f,1)
    delta = 2/(n-1)

    s     = (n,m)
    res   = np.zeros(s)
    # Constructing A
    G1  = sp.diags([1, -2, 1], [-1, 0, 1], (m-2, m-2) )
    Gx  = sp.kron(sp.eye(n-2), G1 )
    G2  = sp.diags([1, -2, 1], [-1, 0, 1], (n-2, n-2) )
    Gy  = sp.kron(G2, sp.eye(m-2) )

    A   = (Gx + Gy)/delta**2
    
    q_vec   = np.ravel(-f[1:-1,1:-1])
    Ain_vec = np.ravel(Ain[1:-1,1:-1])
    
    res_vec = A*Ain_vec - q_vec
    
    res[1:-1,1:-1] = np.reshape(res_vec, (n-2, m-2))
    
    return res



def restrict(Afine): 

#=====================================================
#
#  restriction routine  
#
#=====================================================

    # resolution
    n  = np.size(Afine,0)
    m  = np.size(Afine,1) 
    n2 = int((n+1)/2)
    m2 = int((m+1)/2)

    s  = (n2,m2)
    
    if ((n+1) % 2 != 0) :
        print("n2 not divisible by 2!!!")
    Acoarse = np.zeros(s) 

    # restriction operation
    for i in range(1,n2-1):
        for j in range(1,m2-1):
            Acoarse[i,j] = 1/4*Afine[2*i,2*j]
            Acoarse[i,j] = Acoarse[i,j] + (Afine[2*i-1,2*j] + Afine[2*i+1,2*j] +  Afine[2*i,2*j-1] + Afine[2*i,2*j+1])/8
            Acoarse[i,j] = Acoarse[i,j] + (Afine[2*i-1,2*j-1] + Afine[2*i+1,2*j+1] +  Afine[2*i+1,2*j-1] + Afine[2*i-1,2*j+1])/16
            
    return Acoarse       



def prolong(Acoarse):

#====================================================
#
#  prolongation routine 
#
#====================================================

    # resolution
    n = np.size(Acoarse,0)
    m = np.size(Acoarse,1)
    n2 = 2*n - 1
    m2 = 2*m - 1
    s  = (n2,m2)

    Afine = np.zeros(s)

    #prolongation operation
    for i in range(1,n-1):
        for j in range(1,m-1):
       
            Afine[2*i  ,2*j  ] = Acoarse[i,j]
            Afine[2*i+1,2*j  ] = 1/2*(Acoarse[i,j] + Acoarse[i+1,j])
            Afine[2*i  ,2*j+1] = 1/2*(Acoarse[i,j] + Acoarse[i,j+1])
            Afine[2*i  ,2*j+1] = 1/4*(Acoarse[i,j] + Acoarse[i,j+1] + Acoarse[i+1,j] + Acoarse[i+1,j+1])
            
    return Afine        
            
     