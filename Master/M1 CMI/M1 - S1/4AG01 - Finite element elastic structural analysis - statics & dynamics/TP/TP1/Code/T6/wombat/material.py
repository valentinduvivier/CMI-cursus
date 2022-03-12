# -*- coding: utf-8 -*-
"""
Module for class `Material`

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
import numpy as np

class Material:
    """ Abstract class for material properties"""
    pass

#class LinearElastic(Material):
#    """ Linear elastic material for bars
#    
#    Attributes
#    ----------
#    Young_modulus : float
#        Material Young modulus :math:`E`
#    """
#    def __init__(self,E=1e6):
#        self.Young_modulus = E
#    

class LinearElastic(Material):
    """ Linear elastic material
    
    Attributes
    ----------
    Young_modulus : float
        material Young modulus :math:`E`
    Poisson_coeff : float
        material Poisson coefficient :math:`\\nu` (with :math:`-1<\\nu<1/2`), 
        ignored for :class:`Bar2D <bar2D.Bar2D>`  and :class:`Beam2D <beam2D.Beam2D>`  elements
    rho : float
        material volumetric mass density :math:`\\rho`
    model : {'plane_strain','plane_stress','axi'}
        type of 2D model
    C : ndarray
        elasticity matrix :math:`[C]` shape=(3,3)
    """
    def __init__(self,E=1e6,nu=0.,rho=0.,model="plane_strain"):
        assert (nu<=0.5) and (nu>-1), "Wrong Poisson coefficient"
        self.Young_modulus = E
        self.Poisson_coeff = nu
        self.rho = rho
        self.model = model
        self.C = self.C_matrix()
        
    def compute_lame_coeff(self):
        """Returns Lamé coefficients :math:`\lambda,\mu`"""
        E = self.Young_modulus
        nu = self.Poisson_coeff
        lamb = E*nu/(1+nu)/(1-2*nu)
        mu = E/2./(1+nu)
        return lamb,mu
        
    def from_lame(self,lamb,mu):
        self.Young_modulus = mu*(3*lamb+2*mu)/(lamb+mu)
        self.Poisson_coeff = lamb/2./(lamb+mu)
    
    def C_matrix(self):
        """Compute elasticity matrix :math:`[C]`"""
        if self.model == "plane_strain":
            lamb,mu = self.compute_lame_coeff()
            return np.array([[lamb+2*mu,lamb,0],[lamb,lamb+2*mu,0],[0,0,mu]])
        elif self.model == "plane_stress":
            E = self.Young_modulus
            nu = self.Poisson_coeff
            return E/(1-nu**2)*np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2.]])
    
    def compute_sigzz(self,Eps):
        if self.model == "plane_strain":
            lamb,mu = self.compute_lame_coeff()
            return lamb*(Eps[0]+Eps[1])
        elif self.model == "plane_stress":
            return 0.
