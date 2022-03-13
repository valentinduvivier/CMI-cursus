#!/usr/bin/env python
# coding: utf-8

# # Aspects numériques des modèles à gradient d'endommagement

# In[1]:

from dolfin import *

# get_ipython().run_line_magic('matplotlib', 'inline')
from damage_gradient import *

problem = "homog"
refinement_level = 0
center_thickness = None
bcs_damage = False
unloading = False
export_results = True

## Material properties
#   - Elastic properties
E = 3e3
nu = 0.3
#   - Damage properties
Gc = 3e-3
l0 = Constant(0.1)
sig0 = None
model = "AT1"
# Loading steps
Nincr = 100
# Maximum imposed displacement
Umax = 3e-3

mesh, facets, domains = setup_geometry(problem, refinement_level, center_thickness)

plt.figure()
if problem == "composite":
    pl = plot(domains, cmap="bwr")
    plt.colorbar(pl)
plot(mesh, linewidth=0.5)
plt.show()

# Maximum number of increments in fixed point
Nitermax = 500
# Convergence tolerance for fixed point
tol = 1e-4
mech_params = (E, nu, Gc, l0, sig0, model, Umax, Nincr)
prob_params = (problem, bcs_damage, unloading, export_results, Nitermax, tol)


# In[2]:


results, u, d, sig = solve_problem(mesh, facets, domains, prob_params, mech_params)


# In[ ]:


plot_results(problem, results, E, Gc, l0, model, d, u, sig, mesh, export_results)


# In[ ]:
