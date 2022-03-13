#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 10:23:02 2019

@author: bleyerj
"""
from dolfin import *
from ufl import Max, Min
from mshr import *
import matplotlib.pyplot as plt
import numpy as np

linear_solver = "lu" # or "mumps"

def setup_geometry(problem="homog", refinement_level=0, hole_spacing=0, hole_radius=0.1, aspect_ratio=1):
    L, W, R = 1., 0.5, hole_radius

    if problem != "perforated":
        L, W = 1., 1.
        mesh = RectangleMesh(Point(0, 0), Point(L, W), 2**refinement_level, 1)
    else:
        domain = Rectangle(Point(0, 0), Point(L, W)) \
        - Ellipse(Point(L/2-hole_spacing, 0.), R/aspect_ratio, R, 40) \
        - Ellipse(Point(L/2+hole_spacing, W), R/aspect_ratio, R, 40)
        mesh = generate_mesh(domain, 30)
        # mesh refinement
        for r in range(min(refinement_level, 3)):
            class CentralPart(SubDomain):
                def inside(self, x, on_boundary):
                    return L/2-((1+1/aspect_ratio)*R)/2-hole_spacing <= x[0] <= L/2+((1+1/aspect_ratio)*R)/2+hole_spacing
            to_refine = MeshFunction("bool", mesh, 2)
            CentralPart().mark(to_refine, True)
            mesh = refine(mesh, to_refine)

    # Define boundaries and boundary integration measure
    facets = MeshFunction("size_t", mesh, 1)
    facets.set_all(0)
    class Left(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], 0)
    class Right(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], L)
    Left().mark(facets, 1)
    Right().mark(facets, 2)

    return mesh, facets


def solve_problem(mesh, facets, prob_params, mech_params):
    problem, unloading, export_results, Nitermax, tol = prob_params
    E, nu, R0, alpha, Umax, Nincr = mech_params

    ds = Measure("ds", subdomain_data=facets)

    lmbda = Constant(E*nu/(1+nu)/(1-2*nu))
    mu = Constant(E/2/(1+nu))
    lmbda_ps = 2*mu*lmbda/(lmbda+2*mu) # effective lambda modulus for plane stress
    def eps(v):
        return sym(grad(v))
    def sigma0(v):
        return lmbda_ps*tr(eps(v))*Identity(2) + 2*mu*eps(v)
    def sigma(v, d):
        return (1-d)*sigma0(v)
    def compute_damage(u):
        Y = 0.5*inner(sigma0(u), eps(u))
        return Min((Y/Constant(R0)-1)/Constant(alpha), 1.-1e-8) # min between formula and 1 as d=[0;1]. Moreover, we condiser d=1-eps to make sure in d~=1 the stiffness matrix is still invertible.

    if unloading:
        load_steps = np.concatenate((np.linspace(0, 3*Umax/4, Nincr//2+1),
                                 np.linspace(3*Umax/4, 0, 10)[1:],
                                 np.linspace(0, Umax, Nincr//2+1)[1:]))
    else:
        load_steps = np.linspace(0, Umax, Nincr+1)

    # function space for the displacement
    V_u = VectorFunctionSpace(mesh, "CG", 1)
    # function space for the damage
    V_d = FunctionSpace(mesh, "DG", 0)
    # function space for the stress
    V_sig = TensorFunctionSpace(mesh, "DG", 0)

    Uimp = Expression("t", t=0, degree=0)
    def point(x, on_boundary):
        return near(x[0], 0) and near(x[1], 0)
    bcs = [DirichletBC(V_u.sub(0), Constant(0.), facets, 1),
           DirichletBC(V_u.sub(0), Uimp, facets, 2),
           DirichletBC(V_u, Constant((0, 0)), point, "pointwise")]


    d = Function(V_d, name="Damage")
    dold = Function(V_d, name="Previous damage")
    dlb = Function(V_d, name="Lower bound bound d_n")
    u = Function(V_u, name="Total displacement")
    v = TrialFunction(V_u)
    u_ = TestFunction(V_u)
    sig = Function(V_sig, name="Current stress")

    # bilinear and l form for elasticity problem
    a = inner(eps(v), sigma(u_, dold))*dx
    l = Constant(0)*u_[0]*dx

    Nincr = len(load_steps)-1
    ffile = XDMFFile("results.xdmf")
    ffile.parameters["functions_share_mesh"] = True
    ffile.parameters["flush_output"] = True
    results = np.zeros((Nincr+1, 3))
    t = 0
    for (i, t) in enumerate(load_steps[1:]):
        Uimp.t = t
        print("Increment {:3d}".format(i+1))
        niter = 0
        for niter in range(Nitermax):
            # Solve displacement
            solve(a == l, u, bcs, solver_parameters={"linear_solver":linear_solver})
            # Compute new damage
            d.assign(local_project(Max(dlb, compute_damage(u)), V_d)) # make sure d can only increase by taking max between current d and one of previous step.
            # Compute difference between new and old damage
            nRes = max(d.vector().get_local()-dold.vector().get_local())
            # Update damage
            dold.assign(d)

            print("    Iteration {:3d}: ||Res||={:5e}".format(niter, nRes))
            if nRes < tol:
                break
        else:
            warning("Too many iterations in fixed point algorithm")

        # compute apparent stress and max damage
        results[i+1,:] = (t, 1/assemble(Constant(1)*ds(2, domain=mesh))*assemble(sigma(u, d)[0,0]*ds(2)),
                          max(d.vector().get_local()))

        # update lower bound for irreversibility
        dlb.assign(d)

        # compute stress
        sig.assign(local_project(sigma(u, d), V_sig))

        # Export to Paraview
        if export_results:
            ffile.write(d, t)
            ffile.write(u, t)
            ffile.write(sig, t)

        ffile.close()

    return results, u, d, sig

def plot_results(problem, results, E, R0, alpha, d, u, sig, mesh, export_results=True):
    Ediss = np.trapz(results[:, 1], x=results[:, 0])
    print("Dissipated energy/R0:", Ediss/float(R0))
    sigmax = max(results[:, 1])
    epsmax = results[np.where(results[:, 1]==sigmax)[0][0],0]
    print("Maximum stress:      ", sigmax)
    print("Strain at max stress:", epsmax)

    plt.figure()
    plt.plot(results[:, 0], results[:, 1], '-x', label="FE solution")
    plt.plot([epsmax, epsmax, 0], [0, sigmax, sigmax], "--k", linewidth=1)
    exx = results[:, 0]
    s = E*exx
    exx0 = (2*R0/E)**0.5
    exxd = exx[exx>=exx0]
    dam = 0*exx
    dam[exx>=exx0] = np.minimum(1., (0.5*E*exxd**2-R0)/alpha/R0)
    s[exx>=exx0] = np.maximum(0,(1-dam[exx>=exx0])*E*exxd)
    if problem == "homog":
        plt.plot(exx, s, '--', label="Homogeneous solution")
    plt.xlabel("Uniaxial strain $\epsilon_{xx}$")
    plt.ylabel("Apparent stress $F/S$ (MPa)")
    plt.legend()
    if export_results:
        plt.savefig("stress_strain_curve.pdf")
    plt.show()

    plt.figure()
    plt.plot(results[:, 0], results[:, 2], '-x', label="FE solution")
    if problem == "homog":
        plt.plot(exx, dam, '--', label="Homogeneous solution")
    plt.xlabel("Uniaxial strain $\epsilon_{xx}$")
    plt.ylabel("Maximum damage")
    plt.ylim(0, 1.2)
    plt.legend()
    if export_results:
        plt.savefig("damage_evol_curve.pdf")
    plt.show()

    plt.figure()
    plt.title("Damage field")
    plot(mesh, linewidth=0.1)
    p = plot(d, cmap="bone_r", vmin=0, vmax=1)
    plt.colorbar(p)
    if export_results:
        plt.savefig("damage.png",dpi=300)
    plt.show()

    plt.figure()
    plt.title("Displacement")
    p = plot(u[0])
    plt.colorbar(p)
    if export_results:
        plt.savefig("displacement.png",dpi=300)
    plt.show()

    plt.figure()
    plt.title("$\sigma_{xx}$ stress")
    p = plot(sig[0,0])
    plt.colorbar(p)
    if export_results:
        plt.savefig("axial_stress.png",dpi=300)
    plt.show()


def local_project(v, V, u=None):
    dv = TrialFunction(V)
    v_ = TestFunction(V)
    a_proj = inner(dv, v_)*dx
    b_proj = inner(v, v_)*dx
    solver = LocalSolver(a_proj, b_proj)
    solver.factorize()
    if u is None:
        u = Function(V)
        solver.solve_local_rhs(u)
        return u
    else:
        solver.solve_local_rhs(u)
        return
