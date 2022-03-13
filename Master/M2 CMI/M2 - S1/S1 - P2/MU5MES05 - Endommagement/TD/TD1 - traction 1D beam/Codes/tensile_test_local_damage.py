#!/usr/bin/env python
# coding: utf-8

# # Aspects numériques des modèles d'endommagement locaux
#
#

# In[ ]:


from dolfin import *
from local_damage import *

problem = "homog"
refinement_level = 0
hole_spacing = 0.
hole_radius = 0.2
aspect_ratio = 1
unloading = False
export_results = True

## Material properties
#   - Elastic properties
E = 6e3
nu = 0.3
#   - Damage properties
R0 = 3e-3
alpha = 1
# Loading steps
Nincr = 100
# Maximum imposed displacement
Umax = 3e-3

mesh, facets = setup_geometry(problem, refinement_level, hole_spacing, hole_radius, aspect_ratio)

plt.figure()
plot(mesh, linewidth=0.5)
plt.show()

# Maximum number of increments in fixed point
Nitermax = 200
# Convergence tolerance for fixed point
tol = 1e-4
mech_params = (E, nu, R0, alpha, Umax, Nincr)
prob_params = (problem, unloading, export_results, Nitermax, tol)


# In[ ]:


results, u, d, sig = solve_problem(mesh, facets, prob_params, mech_params)


# In[ ]:


plot_results(problem, results, E, R0, alpha, d, u, sig, mesh, export_results)


# In[ ]:




