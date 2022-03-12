"""
Created on Wed Jun 13 17:39:16 2018
Last modification:
Test Wing 
@author: SophieDartois
@email: sophie.dartois@sorbonne-universite.fr

This file is part of the WomBat finite element code
used for the Civil Engineering Finite Element Course
of Ecole des Ponts ParisTech 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
import matplotlib.pyplot as plt
from wombat import *

# Problem parameters
Yg_mod = 1.
Nu_poi = 0.3
FLoad = 0.0625

# Geometry description and mesh generation
mesh, boundary = call_gmsh('geometry/wing.geo',SolidT3)
right = boundary.get_elem_from_tag("right")
left = boundary.get_elem_from_tag("left")

# Material properties
mat = LinearElastic(E=Yg_mod,nu=Nu_poi,model = "plane_stress")

# Boundary conditions and load
appuis = Connections()
appuis.add_imposed_displ(left,ux=0,uy=0)

forces = ExtForce()
forces.add_distributed_forces(right, fy=FLoad)

# Model construction
model = Model(mesh,mat)

# Assembly phase
K = assembl_stiffness_matrix(model)
L,Ud = assembl_connections(appuis,model)
F = assembl_external_forces(forces,model)

# Solving
U,lamb = solve(K,F,L,Ud)

# Post-processing
list_res = specific_res_treat()

nbFig = 1
coeff = 0.5
analysis_type=mat.model
list_res.treat_stat(mesh)
nbFig = list_res.treat_disp(mesh,coeff,U,appuis,analysis_type,nbFig)
nbFig = list_res.treat_stress(mesh,coeff,U,model,analysis_type,nbFig)

