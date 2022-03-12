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
    Sigeq : Von Mises equivalent stress
    Ep : Potential energy
    
    """    
    
    def __init__(self, UMax=0, VMMax=0):
        self.Max_disp = UMax
        self.Max_VM = VMMax
        self.Sigxx = []
        self.Sigyy = []
        self.Sigxy = []
        self.Sigzz = []
        self.Sigeq = []
        
    def treat_stat(self,mesh):
        print ('')
        print ("nb of nodes in the mesh: ", mesh.nodes.nb_nodes)
        nbelem = len(mesh.elem_list)
        print ("nb of elements in the mesh", nbelem)
        print ('')
        
    def treat_disp(self,mesh,coeff,U,appuis,analysis_type,nbFig):
        
        # Maximum displacement
        Ux = U[0::2]
        Uy = U[1::2]
        normU = (Ux*Ux+Uy*Uy)**(1./2.)
        mU = max(normU)
        print ("Max displacement:"), mU
        
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
            
 
    def treat_stress(self,mesh,coeff,U,model,analysis_type,nbFig):
        
        # Treating and plotting stress fields
        if analysis_type == "plane_strain" or analysis_type == "plane_stress": 
        
            Sigma = stresses(U,model)

            self.Sigxx = Sigma[::4]
            self.Sigyy = Sigma[1::4]
            self.Sigzz = Sigma[2::4]
            self.Sigxy = Sigma[3::4]
                           
            S1 = self.Sigxx-self.Sigyy
            S2 = self.Sigyy-self.Sigzz
            S3 = self.Sigzz-self.Sigxx
            
            # As defined in planar elasticity (plane x,y)
            self.Sigeq = (0.5*(S1**2+S2**2+S3**2+6*self.Sigxy**2))**0.5
            
            store_max = []
            store_max.append(max(self.Sigxx))
            print ("Maximum Sigxx stress:",max(abs(self.Sigxx)))
            store_max.append(max(self.Sigyy))
            print ("Maximum Sigyy stress:",max(abs(self.Sigyy)))
            store_max.append(max(self.Sigxy))
            print ("Maximum Sigxy stress:",max(abs(self.Sigxy)))
            VMmax = max(self.Sigeq)
            print ("Maximum von Mises stress:", VMmax)
            
            for i, j in enumerate(self.Sigeq) :
                if j == VMmax:
                    print ("Max VMstress element :", i)
                    print ("Max VMstress element node coord:",\
                    mesh.elem_list[i].node_coor())
    
            fig = Figure(nbFig,r"$\sigma_{xx}$ stress")
            fig.plot_field(mesh,self.Sigxx)
            fig.show()
    
            fig = Figure((nbFig+1),r"$\sigma_{yy}$ stress")
            fig.plot_field(mesh,self.Sigyy)
            fig.show() 
            
            fig = Figure((nbFig+2),r"$\sigma_{xy}$ stress")
            fig.plot_field(mesh,self.Sigxy)
            fig.show() 
    
            fig = Figure((nbFig+3),"Equivalent von Mises stress")
            fig.plot_field(mesh,self.Sigeq)
            fig.show()
            
            nbFig = nbFig+4
            return nbFig
        
        else: print ("not implemented yet")
    
#Change spatinet TP1
    def treat_energy(self,U,K,F):
                   
            self.Ep = (1./2.)*np.inner(U,K.dot(U))-np.inner(F,U)
            
            print ("Potential energy: ",self.Ep)
#Change spatinet TP1
