from __future__ import division
import numpy as np
import pylab as plt

class solver_WE():
    def __init__(self, x, t, p0, c0, v0=None, rho0=None, go=0):
        """
        class solver_WE : resolution numérique de l'équation des ondes 1D en milieu hétérogène
        Input:
        ------
        x       [np array (Nx,)]           : maillage spatial
        t       [np array (Nt,)]           : maillage temps
        p0      [np array (Nx,)]           : condition initiale en pression
        c0      [np array (Nx,)]           : vitesse de propagation en chaque point du maillage spatial
        v0      [np array (Nx,)-optionnel] : condition initiale en vitesse
        rho0    [np array (Nx,)-optionnel] : maillage spatial
        go      [int]                      : =1 lance la simul; =0 initialise seulement les tableaux

        Output:
        -------
        p       [np array (Nt,Nx)]         : pression pour chaque point espace-temps
        v       [np array (Nt,Nx)]         : vitesse pour chaque point espace-temps
        """
        self.x = x
        self.t = t
        self.p0 = p0
        try:
            if v0==None:
                self.v0 = np.zeros_like(p0)    
        except:
            self.v0=v0
        try:
            if rho0==None:
                self.rho0 = np.ones_like(p0)
        except:
            self.rho0=rho0
        self.c0 = c0
        self.p = np.zeros((len(t),len(x)))
        self.p[0,:]=self.p0
        self.v = np.zeros((len(t),len(x)))
        self.v[0,:]=self.v0
        self.dx = x[1]-x[0]
        self.dt = t[1]-t[0]
        self.cfl = np.max(self.c0) * self.dt/self.dx
        
        # verification de securité
        assert (np.shape(self.v0)==np.shape(self.p0)),"v0 et p0 doivent avoir la meme taille"
        assert (np.shape(self.rho0)==np.shape(self.p0)),"rho0 et p0 doivent avoir la meme taille"
        assert (np.shape(self.c0)==np.shape(self.p0)),"c0 et p0 doivent avoir la meme taille"
        assert (self.cfl<=1), "le nombre cfl=c0.dt/dx doit etre inférieur à 1"
        if go==1:
            self.compute_solution()
        

    def init_solver(self):
        """ initialisation par un schema de Lax-Friedrichs"""
        n = 1
        Nx = np.shape(self.p)[1]
        for j in range(1,Nx-1):
            self.p[n,j] = 0.5 * (self.p[n-1,j+1]+self.p[n-1,j-1]) - self.rho0[j] * self.c0[j]**2 * self.dt / (2 * self.dx) * (self.v[n-1,j+1] - self.v[n-1, j-1])
            self.v[n,j] = 0.5 * (self.v[n-1,j+1]+self.v[n-1,j-1]) - 1/self.rho0[j] * self.dt / (2 * self.dx) * (self.p[n-1,j+1] - self.p[n-1, j-1])

    def step(self):
        """ Leap Frog system """
        Nt = np.shape(self.p)[0]
        Nx = np.shape(self.p)[1]
        J = range(1,Nx-1)
        Jm1 = range(0,Nx-2)
        Jp1 = range(2,Nx)
        for n in range(2,Nt):
            self.p[n,J] = self.p[n-2,J] - self.dt/self.dx * self.rho0[J] * self.c0[J]**2 * (self.v[n-1,Jp1]-self.v[n-1,Jm1])
            self.v[n,J] = self.v[n-2,J] - self.dt/self.dx * 1/self.rho0[J] * (self.p[n-1,Jp1]-self.p[n-1,Jm1])

    def compute_solution(self):
        """ calcule de la solution"""
        self.init_solver()
        self.step()
        
