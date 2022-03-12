#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 18:16:20 2018

Compiles post processing functions

@author: SophieDartois
"""

import numpy as np
import matplotlib.pyplot as plt 

from .post_process import Figure
from .finite_elements_sparse import stresses
import scipy.sparse as sparse

from CalculSigma import *

lin_shape_functions = lambda xi,eta: np.array([xi,eta,1-xi-eta])


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
        self.Sigxx = []
        self.Sigyy = []
        self.Sigxy = []
        self.Sigzz = []
        # create class for Sigma vm
        self.Sigvm = []
        self.Ep = []
        
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
            
 
    def treat_stress(self,mesh,coeff,U,model,analysis_type,nbFig):
        
        # Treating and plotting stress fields
        if analysis_type == "plane_strain" or analysis_type == "plane_stress": 
        
            Sigma = stresses(U,model)
            
            self.Sigxx = Sigma[::4]
            self.Sigyy = Sigma[1::4]
            self.Sigzz = Sigma[2::4]
            self.Sigxy = Sigma[3::4]
            
            # Sigma vm
            self.Sigvm=(1/np.sqrt(2))*np.sqrt(((self.Sigxx-self.Sigyy)**2 + (self.Sigxx-self.Sigzz)**2 +(self.Sigyy-self.Sigzz)**2  +6*(self.Sigxy)**2 ))
            
            store_max = []
            store_max.append(max(self.Sigxx))
            print ("Maximum Sigxx stress: "+self.fmt.format(max(abs(self.Sigxx))))
            store_max.append(max(self.Sigyy))
            print ("Maximum Sigyy stress: "+self.fmt.format(max(abs(self.Sigyy))))
            store_max.append(max(self.Sigxy))
            print ("Maximum Sigxy stress: "+self.fmt.format(max(abs(self.Sigxy))))
            
            store_max.append(max(self.Sigvm))
            print ("Maximum Sigvm stress: "+self.fmt.format(max(abs(self.Sigvm))))
            
    
            fig = Figure(nbFig,r"$\sigma_{xx}$ stress")
            fig.plot_field(mesh, self.Sigxx)
            fig.show()
    
            fig = Figure((nbFig+1),r"$\sigma_{yy}$ stress")
            fig.plot_field(mesh, self.Sigyy)
            fig.show() 
            
            fig = Figure((nbFig+2),r"$\sigma_{xy}$ stress")
            fig.plot_field(mesh, self.Sigxy)
            fig.show() 
            
            # plot Sigma vm
            fig = Figure((nbFig+3),r"$\sigma_{vm}$ stress")
            fig.plot_field(mesh, self.Sigvm)
            fig.show() 

            nbFig = nbFig+4
            return nbFig
        
        else: print ("not implemented yet")
    
    def treat_potential_energy(self,U,K,F):
        "definition of the potential energy"
        
        W_el = np.dot(U,K.dot(U))/2.
        Phi_ext = np.dot(F,U)
        Ep = W_el-Phi_ext
        
        print("Potential energy : "+self.fmt.format(Ep))

                
        
    def extract_Sig(self, mesh, V, theta):
        
        """ Sigma extraction """
        
        nodes_rad_coord = []
        V_extract       = []
        e0 = mesh.elem_list[0]
        
        try: 
            ngauss = e0.ngauss
        except:
            ngauss = 1
            
        if V.shape[0] == mesh.Nel:      # case of constant field/elem
            print('constant field')    
            
            for(j,e) in enumerate (mesh.elem_list):
                node_coords = e.node_coor()
                    
                for ind_co in range (node_coords.shape[0]):
                    co = node_coords[ind_co]
                        
                    if np.abs(np.arctan2(co[1],co[0])-theta)<0.02:
                        co_rad=np.sqrt(co[0]**2 + co[1]**2)
                        
                        nodes_rad_coord.append(co_rad)
                        V_extract.append(V[j])
            
            
        elif V.shape[0] == ngauss*mesh.Nel and ngauss > 1:  # case of non uniform field/elem
            print ("Case of non uniform field")
            
            for (j,e) in enumerate(mesh.elem_list):
                node_coords = e.node_coor()[:3]
                N = np.array([lin_shape_functions(e.ag[i,0],e.ag[i,1]) for i in range(ngauss)])
                V_node = np.linalg.lstsq(N,V[ngauss*j:ngauss*(j+1)])[0]
                
                for ind_co in range (node_coords.shape[0]):
                        co = node_coords[ind_co]
                        
                        if np.abs(np.arctan2(co[1],co[0])-theta ) < 0.02:
                            co_rad=np.sqrt(co[0]**2+co[1]**2)
                            
                            nodes_rad_coord.append(co_rad)
                            V_extract.append(V_node[ind_co])
            
        else: # case of nodal field for T3 or constant field/element
            print("case of nodal field for T3 or constant field/element")
        
        return (nodes_rad_coord, V_extract)
    
    def plot_extract_Sig(self, mesh, V_num, R, V_ana, r, theta, idx, idx2):

        
        plt.figure()
        plt.plot(R, V_num, 'r.', label='Numerical stress')
        plt.plot(r, V_ana, 'b-', label='Analytical stress')
        
        if idx == 1:
            plt.title(r'$\sigma_{rr}$ at $\theta$ = %.2f°' %(theta))
        else:
            plt.title(r'$\sigma_{yy}$ at $\theta$ = %.2f°' %(theta))
        plt.xlabel('r [m]')
        plt.ylabel(f'Stress [Pa]')
        
        plt.grid('True')
        plt.legend()

        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        
        # -----------------------------------------------------------------
        
