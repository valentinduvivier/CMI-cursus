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

import numpy as np
import matplotlib.pyplot as plt

from wombat import *
from CalculSigma import *

# Material parameters
Yg_mod = 135*10**9      # [Pa] == 35 GPa
Nu_poi = 0.35           # [1]
rho_c  = 2300           # [kg.m-3]

g      = 9.81            # [m.s-2]

# Geometry parameters
Length = 40*10**-2/2      # [m]
b      = 2.5*10**-2     # [m]
a      = 2*b 

# Forces
fs     = 10**6             # [N]
Weight = - 1.*rho_c*g   # [N]

# Geometry description and mesh generation
mesh, boundary = call_gmsh('geometry/PlateHole_Symmetry.geo', SolidT3)

left   = boundary.get_elem_from_tag("left")
right  = boundary.get_elem_from_tag("right")

top    = boundary.get_elem_from_tag("top")
bottom = boundary.get_elem_from_tag("bottom")

hole   = boundary.get_elem_from_tag("hole")
plate  = boundary.get_elem_from_tag("plate_ellipse")

# Material properties
mat = LinearElastic(E=Yg_mod, nu=Nu_poi, rho=rho_c, model = "plane_stress")

# Boundary conditions and load
appuis = Connections()
appuis.add_imposed_displ(bottom, uy=0)
appuis.add_imposed_displ(left, ux=0)

forces = ExtForce()
forces.add_distributed_forces(top, fy=fs)

# no force on bottom as working on symmetry
#forces.add_distributed_forces(bottom, fy=-F_s)

forces.add_volume_forces(plate, fy=Weight)

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

# ------------------------------------------------------------------------------

# Analysis sigma (comparison numeric/analytic)

r = np.linspace(a, Length/1.5, 55)

theta = 0#np.pi/2
idx = 2
idx2 = 'r'
(R, V_num) = list_res.extract_Sig(mesh, list_res.Sigyy, theta)
V_ana = sig_tt(fs, a, r, theta)
list_res.plot_extract_Sig(mesh, V_num, R, V_ana, r, theta, idx, idx2)

# ------------------------------------------------------------------------------

# K_G

N     = 50
r     = 0.025
theta = np.linspace(0, np.pi/2, N)

Stresses = np.zeros((3, N))
Siggg = [list_res.Sigxx, list_res.Sigxy, list_res.Sigyy]

sigma_G = fs
sigma_M = np.max(Siggg)
K_G     = sigma_M/sigma_G

print(K_G)
