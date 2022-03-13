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

linear_solver = "lu"  # or "mumps"


def setup_geometry(
    problem="homog",
    refinement_level=0,
    thickness=0.1,
):
    if problem == "homog":
        L, W = 1.0, 0.1
        if refinement_level > 0:
            mesh = generate_mesh(
                Rectangle(Point(0, 0), Point(L, W)), 2 ** (refinement_level)
            )
        else:
            mesh = RectangleMesh(Point(0, 0), Point(L, W), 2 ** refinement_level, 1)
    elif problem == "composite":
        L, W = 1, 0.2
        domain = Rectangle(Point(0, 0), Point(L, W))
        domain.set_subdomain(
            1,
            Rectangle(Point(0, W / 2 - thickness / 2), Point(L, W / 2 + thickness / 2)),
        )
        mesh = generate_mesh(domain, 60 + 30 * (refinement_level))
    else:
        raise ValueError("Wrong problem name '{}' is not valid.".format(problem))

    # Define boundaries and boundary integration measure
    facets = MeshFunction("size_t", mesh, 1)
    facets.set_all(0)

    class Left(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], 0)

    class Right(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], L)

    class Bottom(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], 0)

    class Top(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], W)

    Left().mark(facets, 1)
    Right().mark(facets, 2)
    Bottom().mark(facets, 3)
    Top().mark(facets, 4)

    domains = MeshFunction("size_t", mesh, 2, mesh.domains())

    return mesh, facets, domains


def cell_project(values_dict, V, dx):
    """Used to define cell-variable coefficients."""
    v_, dv = TestFunction(V), TrialFunction(V)
    a_form = inner(v_, dv) * dx
    L_form = sum([inner(v_, val) * dx(reg) for (reg, val) in values_dict.items()])
    u = Function(V)
    A = assemble(a_form)
    b = assemble(L_form)
    solve(A, u.vector(), b)
    return u


def solve_problem(mesh, facets, domains, prob_params, mech_params):
    problem, bcs_damage, unloading, export_results, Nitermax, tol = prob_params
    E, nu, Gc, l0, sig0, model, Umax, Nincr = mech_params

    ds = Measure("ds", subdomain_data=facets)
    dx = Measure("dx", subdomain_data=domains)

    lmbda = Constant(E * nu / (1 + nu) / (1 - 2 * nu))
    mu = Constant(E / 2 / (1 + nu))
    lmbda_ps = (
        2 * mu * lmbda / (lmbda + 2 * mu)
    )  # effective lambda modulus for plane stress
    kres = Constant(1e-6)  # residual stiffness

    # function space for the displacement
    V_u = VectorFunctionSpace(mesh, "CG", 1)
    # function space for the damage
    V_d = FunctionSpace(mesh, "CG", 1)
    # function space for the stress
    V_sig = TensorFunctionSpace(mesh, "DG", 0)
    # function space for material properties
    V0 = FunctionSpace(mesh, "DG", 0)

    if problem == "homog":
        Gc_h = Constant(Gc)
    elif problem == "composite":
        Gc_h = cell_project({0: 1e6 * Gc, 1: Gc}, V0, dx)
        # we perturb the values by a Gaussian distribution
        np.random.seed(15012017)
        Gc_h.vector()[:] *= np.random.normal(1, 0.2, (V0.dim(),))

    def degrad(d):
        if model in ["AT1", "AT2"]:
            return (1 - d) ** 2 + kres
        elif model == "DM":
            lc = E * Gc / sig0**2
            K  = 2*lc/(l0*cw) - 1
            a  = (1 - d) ** 2 / (1 + K*w(d))
            return a + kres
        else:
            raise ValueError("Wrong model name:", model)

    def eps(v):
        return sym(grad(v))

    def sigma0(v):
        return lmbda_ps * tr(eps(v)) * Identity(2) + 2 * mu * eps(v)

    def sigma(v, d):
        return degrad(d) * sigma0(v)

    if unloading:
        load_steps = np.concatenate(
            (
                np.linspace(0, 3 * Umax / 4, Nincr // 2 + 1),
                np.linspace(3 * Umax / 4, 0, 20)[1:],
                np.linspace(0, Umax, Nincr // 2 + 1)[1:],
            )
        )
    else:
        load_steps = np.linspace(0, Umax, Nincr + 1)

    Uimp = Expression("t", t=0, degree=0)

    def point(x, on_boundary):
        return near(x[0], 0) and near(x[1], 0)

    bcs = [
        DirichletBC(V_u.sub(0), Constant(0.0), facets, 1),
        DirichletBC(V_u.sub(0), Uimp, facets, 2),
        DirichletBC(V_u, Constant((0, 0)), point, "pointwise"),
    ]

    u = Function(V_u, name="Total displacement")
    v = TrialFunction(V_u)
    u_ = TestFunction(V_u)
    sig = Function(V_sig, name="Current stress")

    d = Function(V_d, name="Damage")
    dold = Function(V_d, name="Previous damage")
    dub = Function(V_d, name="Upper bound d=1")
    dub.interpolate(Constant(1.0))
    dlb = Function(V_d, name="Lower bound bound d_n")
    d_ = TestFunction(V_d)
    dd = TrialFunction(V_d)
    # boundary conditions for damage problem
    if bcs_damage:
        bc_d = [
            DirichletBC(V_d, Constant(0.0), facets, 1),
            DirichletBC(V_d, Constant(0.0), facets, 2),
        ]
    else:
        bc_d = []
    for bc in bc_d:
        bc.apply(dub.vector())

    if model == "AT1":
        cw = Constant(8 / 3.0)
        w  = lambda d: d
    elif model == "AT2":
        cw = Constant(2.0)
        w  = lambda d: d ** 2
    elif model == "DM":
        cw = Constant(np.pi)
        w  = lambda d: 1 - (1 - d) ** 2

    # bilinear and l form for elasticity problem
    a = inner(eps(v), sigma(u_, dold)) * dx
    l = Constant(0) * u_[0] * dx

    # total energy for the damage problem
    psi = 0.5 * inner(sigma(u, d), eps(u))
    elastic_energy = psi * dx

    fracture_energy = Gc_h / cw * (w(d) / l0 + l0 * dot(grad(d), grad(d))) * dx
    total_energy = elastic_energy + fracture_energy
    # first derivative of energy with respect to d
    F_dam = derivative(total_energy, d, d_)
    # second derivative of energy with respect to d
    J_dam = derivative(F_dam, d, dd)
    # Definition of the optimisation problem with respect to d
    class DamageProblem(OptimisationProblem):
        def f(self, x):  # computation of objective function
            d.vector()[:] = x
            return assemble(total_energy)

        def F(self, b, x):  # computation of first derivative
            d.vector()[:] = x
            assemble(F_dam, b)

        def J(self, A, x):  # computation of second derivative
            d.vector()[:] = x
            assemble(J_dam, A)

    # damage optimization solver
    solver = PETScTAOSolver()
    solver.parameters["method"] = "tron"
    solver.parameters["line_search"] = "gpcg"
    solver.parameters["maximum_iterations"] = 500
    solver.parameters["linear_solver"] = linear_solver
    solver.parameters["gradient_absolute_tol"] = tol
    solver.parameters["gradient_relative_tol"] = tol
    solver.parameters["gradient_t_tol"] = tol

    Nincr = len(load_steps) - 1
    ffile = XDMFFile("results.xdmf")
    ffile.parameters["functions_share_mesh"] = True
    ffile.parameters["flush_output"] = True
    results = np.zeros((Nincr + 1, 4))
    t = 0
    for (i, t) in enumerate(load_steps[1:]):
        Uimp.t = t
        print("Increment {:3d}".format(i + 1))
        niter = 0
        for niter in range(Nitermax):
            # Solve displacement
            solve(a == l, u, bcs, solver_parameters={"linear_solver": linear_solver})
            # Compute new damage
            solver.solve(DamageProblem(), d.vector(), dlb.vector(), dub.vector())
            # Compute difference between new and old damage
            nRes = max(d.vector().get_local() - dold.vector().get_local())
            # Update damage
            dold.assign(d)

            print("    Iteration {:3d}: ||Res||={:5e}".format(niter, nRes))
            if nRes < tol:
                break
        else:
            warning("Too many iterations in fixed point algorithm")

        # Update lower bound to account for irreversibility
        dlb.assign(d)

        resultant = (
            1
            / assemble(Constant(1) * ds(2, domain=mesh))
            * assemble(sigma(u, d)[0, 0] * ds(2))
        )
        results[i + 1, :] = (
            t,
            resultant,
            assemble(elastic_energy) / float(Gc),
            assemble(fracture_energy) / float(Gc),
        )

        sig.assign(local_project(sigma(u, d), V_sig))

        # Export to Paraview
        if export_results:
            ffile.write(d, t)
            ffile.write(u, t)
            ffile.write(sig, t)

            np.savetxt(
                "results_data.csv",
                np.asarray(results),
                header="Uimp, Apparent stress, Elastic energy/Gc, Fracture energy/Gc",
                delimiter=",",
            )

        ffile.close()

    return results, u, d, sig


def plot_results(
    problem, results, E, Gc, l, model, d, u, sig, mesh, export_results=True
):
    print("Dissipated energy/Gc:", results[-1, 3])
    sigmax = max(results[:, 1])
    epsmax = results[np.where(results[:, 1] == sigmax)[0][0], 0]
    print("Maximum stress:      ", sigmax)
    print("Strain at max stress:", epsmax)

    if problem != "shear":
        xs = np.linspace(0, 1.0, 1000)
        plt.figure()
        plt.plot(
            xs,
            [d(x, max(mesh.coordinates()[:, 1]) / 2) for x in xs],
            "-",
            label="damage profile",
        )
        plt.xlabel("x")
        plt.ylabel("Damage")
        plt.ylim(0, 1.0)
        plt.title("Damage profile along $y=W/2$")
        if export_results:
            plt.savefig("damage_profile.pdf")

    plt.figure()
    plt.plot(results[:, 0], results[:, 1], "-x", label="FE solution")
    plt.plot([epsmax, epsmax, 0], [0, sigmax, sigmax], "--k", linewidth=1)
    exx = results[:, 0]
    dam = 0 * exx
    s = E * exx
    if model == "AT1":
        cw = 8 / 3.0
        exx0 = (float(Gc) / float(l) / cw / float(E)) ** 0.5
        exxd = exx[exx >= exx0]
        dam[exx >= exx0] = 1 - (exx0 / exxd) ** 2
        s[exx >= exx0] = (1 - dam[exx >= exx0]) ** 2 * E * exxd
    elif model == "AT2":
        cw = 2.0
        exx0 = (float(Gc) / float(l) / cw / float(E)) ** 0.5
        dam[1:] = 1 / (1 + 2 * exx0 ** 2 / exx[1:] ** 2)
        s = (1 - dam) ** 2 * E * exx
    if problem == "homog" and model in ["AT1", "AT2"]:
        plt.plot(exx, s, "--", label="Homogeneous solution")
    plt.xlabel("Uniaxial strain $\epsilon_{xx}$")
    plt.ylabel("Apparent stress $F/S$ (MPa)")
    plt.legend()
    if export_results:
        plt.savefig("stress_strain_curve.pdf")
    plt.show()

    plt.figure()
    plt.plot(results[:, 0], results[:, 3], "-x", label="fracture (FE solution)")
    plt.plot(results[:, 0], results[:, 2], "-x", label="elastic (FE solution)")
    plt.plot(
        results[:, 0], results[:, 2] + results[:, 3], "-x", label="total (FE solution)"
    )

    plt.xlabel("Uniaxial strain $\epsilon_{xx}$")
    plt.ylabel("Energies/Gc")
    plt.legend()
    if export_results:
        plt.savefig("fracture_energy_evol_curve.pdf")
    plt.show()

    plt.figure()
    plt.title("Damage field")
    plot(mesh, linewidth=0.1)
    p = plot(d, cmap="bone_r", vmin=0, vmax=1)
    plt.colorbar(p)
    if export_results:
        plt.savefig("damage.png", dpi=300)
    plt.show()

    plt.figure()
    plt.title("Displacement")
    p = plot(u[0])
    plt.colorbar(p)
    if export_results:
        plt.savefig("displacement.png", dpi=300)
    plt.show()


def local_project(v, V, u=None):
    dv = TrialFunction(V)
    v_ = TestFunction(V)
    a_proj = inner(dv, v_) * dx
    b_proj = inner(v, v_) * dx
    solver = LocalSolver(a_proj, b_proj)
    solver.factorize()
    if u is None:
        u = Function(V)
        solver.solve_local_rhs(u)
        return u
    else:
        solver.solve_local_rhs(u)
        return
