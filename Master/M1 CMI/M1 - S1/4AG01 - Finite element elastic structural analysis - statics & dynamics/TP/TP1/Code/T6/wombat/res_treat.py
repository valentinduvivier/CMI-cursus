#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 18:16:20 2018

Compiles post processing functions

@author: SophieDartois
"""

import numpy as np
#from post_process import *
from .post_process import Figure
from .finite_elements_sparse import stresses
import scipy.sparse as sparse


class specific_res_treat:
    
    """ Specific results treatment. after solving. 
    To be enriched with other functions and  processes
    
    Attributes
    ----------
    Max_disp : float
        Maximum norm of the displacement vector field.
    Max_VM : float
        Maximum equivalent Von Mises stress
    Sigxx, Sigyy, Sigxy, Sigzz : Components of the stress tensor
    Ep : Potential energy
    """
    
    def __init__(self, UMax=0):
        self.fmt = "{:.2e}" # format to dispay numerical data : 2-digits scientific notation
        self.Max_disp = UMax
        
        self.Sigxx  = []
        self.Sigyy  = []
        self.Sigxy  = []
        self.Sigzz  = []
        
        self.Sig_VM = []
        
        # Energy term
        self.Ep     = []
        
    def treat_stat(self,mesh):
        print ('')
        print ("nb of nodes in the mesh: {:d}".format(mesh.nodes.nb_nodes))
        nbelem = len(mesh.elem_list)
        print ("nb of elements in the mesh: {:d}".format(nbelem))
        print ('')
        
    def treat_disp(self,mesh,coeff,U,appuis,analysis_type,nbFig):
        
        # Maximum displacement
        Ux = U[0::2]
        Uy = U[1::2]
        normU = (Ux*Ux+Uy*Uy)**(1./2.)
        mU = max(normU)
        print ("Max displacement: "+self.fmt.format(mU))
        
        for i, j in enumerate(normU) :
            if j == mU:
                print ("Max displacement node:", i)
                print ("Max displacement node coord:", mesh.nodes.coor[i])
        
        # Displaying deformed shape

        if analysis_type == "plane_strain" or analysis_type == "plane_stress": 
        
            fig = Figure(nbFig,"Deformed shape") 
            fig.plot_def(mesh,U,coeff)
            fig.plot_bc(mesh,appuis)
            fig.show()
            
            nbFig = nbFig+1
            return nbFig
        
        else: print ("not implemented yet")
            
 
    def treat_stress(self, mesh, coeff, U, model, analysis_type, nbFig, E, nu):
        
        # Treating and plotting stress fields
        if analysis_type == "plane_strain" or analysis_type == "plane_stress": 
        
            Sigma = stresses(U, model, E, nu)

            self.Sigxx = Sigma[::4]
            self.Sigyy = Sigma[1::4]
            self.Sigzz = Sigma[2::4]
            self.Sigxy = Sigma[3::4]
            
            self.Sig_VM = (1/np.sqrt(2))*np.sqrt((self.Sigxx - self.Sigyy)**2 + (self.Sigxx - self.Sigzz)**2 + (self.Sigyy - self.Sigzz)**2 + 6*(self.Sigxy)**2)
            
            store_max = []
            store_max.append(max(self.Sigxx))
            print ("Maximum Sigxx stress: "+self.fmt.format(max(abs(self.Sigxx))))
            store_max.append(max(self.Sigyy))
            print ("Maximum Sigyy stress: "+self.fmt.format(max(abs(self.Sigyy))))
            store_max.append(max(self.Sigxy))
            print ("Maximum Sigxy stress: "+self.fmt.format(max(abs(self.Sigxy))))
            store_max.append(max(self.Sig_VM))
            print ("Maximum SigVM stress: "+self.fmt.format(max(abs(self.Sig_VM))))
            
            print("\nSigVM stress vector")
            idx = np.where(self.Sig_VM==max(self.Sig_VM))
            print(self.Sig_VM[idx])

            fig = Figure(nbFig,r"$\sigma_{xx}$ stress")
            fig.plot_field(mesh,self.Sigxx)
            fig.show()
    
            fig = Figure((nbFig+1),r"$\sigma_{yy}$ stress")
            fig.plot_field(mesh,self.Sigyy)
            fig.show() 
            
            fig = Figure((nbFig+2),r"$\sigma_{xy}$ stress")
            fig.plot_field(mesh,self.Sigxy)
            fig.show() 
    
            fig = Figure((nbFig+3),r"$\sigma_{VM}$ stress (mean stress)")
            fig.plot_field(mesh,self.Sig_VM)
            fig.show() 
            
            nbFig = nbFig+4
            return nbFig
        
        else: print ("not implemented yet")
    
    def treat_potential_energy(self, F, K, U):
    
        """ Assembly procedure of the external forces vector
    
            Parameters
            ----------
            forces
                :class:`ExtForce <forces.ExtForce>` object
            model
                :class:`Model <model.Model>` object
    
            Returns
            -------
            F : ndarray
                :math:`\{F\}` vector of equivalent external forces shape=(Nd,)
                where Nd number of dofs
        """
        
        Ep = .5*np.dot(U,K.dot(U)) - np.dot(F, U)
        
        print("\nEp - potential energy")
        print(f"{Ep:.2f} J")
    