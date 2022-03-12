"""
Created on Thu Dec 03 14:59 2020
Last modification:
Test Dam 
@author: DuvivierValentin
@email: valentin.duvivier@etu.sorbonne-universite.fr

This file is an application of WomBat finite element code
used for Gravity dam analysis through T3 element in a first
part and then through a T6 element

"""

import matplotlib.pyplot as plt
from wombat import *

# Problem parameters
Yg_mod = 35*10**9   # [Pa] == 35 GPa
Nu_poi = 0.2        # [1]
rho_c  = 2300       # [kg.m-3]

# Data FLoad - case 1 : homogeneous pressure
patm  = 10**5  # [Pa]
rho_w = 10**3  # [kg.m-3]
g     = 9.81   # [m.s-2]

h     = 15     # [m]

FLoad  = .5*rho_w*g*h     # [N]
Weight = - 1.*rho_c*g

# Geometry description and mesh generation
mesh, boundary = call_gmsh('geometry/GeometryDam_T3.geo', SolidT3)

left   = boundary.get_elem_from_tag("left")
bottom = boundary.get_elem_from_tag("bottom")
dam    = mesh.get_elem_from_tag("dam")

# Material properties
mat = LinearElastic(E=Yg_mod,nu=Nu_poi,rho=rho_c,model = "plane_stress")

# Boundary conditions and load
appuis = Connections()
appuis.add_imposed_displ(bottom,ux=0,uy=0)

forces = ExtForce()
forces.add_distributed_forces(left, fx=FLoad)
forces.add_volume_forces(dam, fy=Weight)

# Model construction
model = Model(mesh,mat)

# Assembly phase
K = assembl_stiffness_matrix(model)
L,Ud = assembl_connections(appuis, model)
F = assembl_external_forces(forces, model)

# Solving
U, lamb = solve(K, F, L, Ud)

# Post-processing
list_res = specific_res_treat()

nbFig = 1
coeff = 0.5
analysis_type=mat.model
list_res.treat_stat(mesh)
nbFig = list_res.treat_disp(mesh, coeff, U, appuis, analysis_type, nbFig)
nbFig = list_res.treat_stress(mesh, coeff, U, model, analysis_type, nbFig)

# Energy calculation
list_res.treat_potential_energy(F, K, U)

import numpy as np
Ep = .5*np.dot(U,K.dot(U)) - np.dot(F, U)
print([mesh.nodes.nb_nodes, Ep])

# ------------------------------------------------

Sigma = stresses(U, model)

Sigxx = Sigma[::4]
Sigyy = Sigma[1::4]
Sigzz = Sigma[2::4]
Sigxy = Sigma[3::4]
            
Sig_VM = (1/np.sqrt(2))*np.sqrt((Sigxx - Sigyy)**2 + (Sigxx - Sigzz)**2 + (Sigyy - Sigzz)**2 + 6*(Sigxy)**2)
    
print([mesh.nodes.nb_nodes, np.max(Sig_VM)])