# -*- coding: utf-8 -*-
"""
Module for class `Beam2D`

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
from .generic_element import *

class Beam2D(Segment):
    """ A 2D Euler-Bernoulli beam element
    
        Beam2D is a Segment-type element (2 end nodes) in 2D
        with 3 degrees of freedom/node :\n 
        
        - **Kinematics**: horizontal, vertical displacement \
        and rotation with a cubic interpolation inside the element  \
        and with Navier-Bernoulli hypothesis :math:`\\theta_z = \\dfrac{du_y}{dx}`
        
        .. math:: \{U\}=\\langle u_x^1,u_y^1,\\theta_z^1,u_x^2,u_y^2,\\theta_z^2\\rangle^T 
            
        - **Strains**: 
            + axial strain :math:`\epsilon`  (constant)
            + bending curvature :math:`\chi` (linear)
        - **Stresses**: 
            + normal force :math:`N` (constant)
            + bending moment :math:`M` (linear)
            + shear force :math:`V` (constant)
    """
    def __init__(self,node_list,tag=1):
        """
        Parameters
        ----------
        
        node_list : list
            list containing two nodes
        tag : int,str
            tag of physical group
        """
        Segment.__init__(self,node_list,tag)
        self.el_dof = 6
        self.node_dof = 3
        self.nb_stresses = 4
        self.kin_field_names = ['U_x','U_y','Theta_z']
        self.strain_field_names = ['eps','chi']
        self.int_forces_field_names = ['N','M1','M2','V']
        self.ext_forces_field_names = ['F_x','F_y','C_z']
    
    def rotation_matrix(self):
        """
        Rotation matrix :math:`[R]` from global to local frame 
        
        shape = (6,6)
        """
        T = self.nodes.coor
        tang = T[1,:]-T[0,:]
        L = self.measure()
        t = tang/L
        
        # rotation matrix [r]{U_x,U_y,Theta_z} = {U_t,U_n,Theta_z}
        r = np.array([[t[0],t[1],0],[-t[1],t[0],0],[0,0,1]])
        # R is such that [R]{U_x1,U_y1,U_x2,U_y2} = {U_t1,U_t2}
        return np.kron(np.eye(2),r) # get a block matrix by repeating [r] twice along the diagonal, 
    
    def shape_functions(self,x):
        """ Cubic shape functions of the Hermite beam element

        Parameters
        ----------
        x : float
            position along the beam axis :math:`x\in[0;L]`
        
        Returns
        -------
        N : ndarray shape=(4,)
            shape functions :math:`[N_i(x)]` evaluated at `x`
        DN : ndarray shape=(4,)
            shape functions first derivatives :math:`[N'_i(x)]` evaluated at `x`
        D2N : ndarray shape=(4,)
            shape functions second derivatives :math:`[N''_i(x)]` evaluated at `x`
            
        """
        L = self.measure()
        N = np.array([1-3*x**2/L**2+2*x**3/L**3,
                      x-2*x**2/L+x**3/L**2,
                      3*x**2/L**2-2*x**3/L**3,
                      -x**2/L+x**3/L**2])
        DN = np.array([6*x/L**2+6*x**2/L**3, 
                       1-4*x/L+3*x**2/L**2, 
                       6*x/L**2-6*x**2/L**3, 
                       -2*x/L+3*x**2/L**2])
        D2N = np.array([6*(-1 + 2*x/L)/L**2,
                        2*(-2 + 3*x/L)/L, 
                        6*(1 - 2*x/L)/L**2,
                        2*(-1 + 3*x/L)/L])
        return N,DN,D2N
    
    def compute_Be_matrix(self,x):
        """ Strain :math:`[B]` matrix such that
        
        .. math:: [B(x)]\{U_e\} = \\begin{Bmatrix} \delta \\\\ \chi(x) \end{Bmatrix}
        
        Parameters
        ----------
        x : float
            position along the beam axis :math:`x\in[0;L]`
            
        """
        L = self.measure()
        N,DN,D2N = self.shape_functions(x)
        B = np.zeros((2,6))
        B[0,[0,3]] = np.array([-1,1])/L
        B[1,[1,4]] = D2N[::2]
        B[1,[2,5]] = D2N[1::2]
        return B
        
        
    def elementary_stiffness(self,mat,sect):
        """ Elementary stiffness matrix :math:`[K_e]` in global frame shape=(12,12)""" 
        L = self.measure()
        R = self.rotation_matrix()
        
        E = mat.Young_modulus
        S = sect.area
        I = sect.inertia
        
        # elementary normal force stiffness matrix in local frame [Ke_normal]{U_t1,U_t2} = {N1,N2}
        Ke_normal_loc = E*S/L*(np.array([[1,-1],[-1,1]]))
        # elementary bending stiffness matrix in local frame [Ke_bend] = {U_n1,Theta_z2,U_n1,Theta_y2}={V1,Mz_1,V2,Mz_2}
        Ke_bend_loc = E*I/L**3*(np.array([[12,6*L,-12,6*L], 
                                         [6*L,4*L**2,-6*L,2*L**2],
                                         [-12,-6*L,12,-6*L],
                                         [6*L,2*L**2,-6*L,4*L**2]]))
                                                        
        Ke_loc = np.zeros((6,6))
        # add normal force stiffness to at rows/columns corresponging to U_t1,U_t2
        Ke_loc[np.ix_([0,3],[0,3])] = Ke_normal_loc
        # add bending moment stiffness to at rows/columns corresponging to U_n1,Theta_z1,U_n2,Thtea_z2
        Ke_loc[np.ix_([1,2,4,5],[1,2,4,5])] = Ke_bend_loc
        
        # elementary stiffness matrix in global frame [Ke_glob] = [R]^T*[Ke_loc]*[R]
        Ke_glob = np.dot(np.dot(R.T,Ke_loc),R)
                   
        return Ke_glob
        
    def elementary_thermal_vector(self,mat,sect,dilat):        
        """ Elementary force vector induced by a thermal strain
        
        Parameters
        ----------
        dilat : float
            uniform thermal dilatation :math:`\\delta_{ther}` inside the element
        """
        T = self.nodes.coor
        tang = T[1,:]-T[0,:]
        L = self.measure()
        tang = tang/L
        
        E = mat.Young_modulus
        S = sect.area
        fther_el = np.zeros((self.el_dof,))
        fther_el[[0,1,3,4]] = E*S*dilat*np.hstack((-tang,tang))
        return fther_el
    
    
    def elementary_distributed_forces(self,el_force):
        """ Elementary force vector for uniform distributed loading
        
        Parameters
        ----------
        el_force = [fx,fy,cz] : array,list
            contains uniformly distributed forces and couple :math:`(f_x,f_y,c_z)`           
        """
        fx,fy,cz = el_force
            
        L = self.measure()
        R = self.rotation_matrix()
        r = R[:3,:3]
        
        P = np.zeros((3,6))
        P[0,[0,3]] = np.array([L/2.,L/2.])
        P[1,[1,2,4,5]] = np.array([L/2.,L**2/12.,L/2.,-L**2/12.])
        P[2,[1,2,4,5]] = np.array([-1,0,1,0])
        
        return np.dot(np.dot(np.dot(np.array([fx,fy,cz]),r.T),P),R)
        
    def deformation(self,Ue,m=21):
        """ Interpolation of the deformed element
        
        Parameters
        ----------
        Ue : ndarray
            nodal displacement of the current elements 
        m : int
            number of points used to interpolate the deformed configurations
        
        Returns
        -------
        x_def,y_def : ndarray
            returns deformed position of m points along the element
        """
        s = np.linspace(-1,1,m)
        x = self.node_coor()[:,0]
        y = self.node_coor()[:,1]
        tang = np.array([x[1]-x[0],y[1]-y[0]])
        L = self.measure()
        t = tang/L
        R = self.rotation_matrix()
        
        Uloc = np.dot(R,Ue)
        Ut = Uloc[0::self.node_dof]
        Un = Uloc[1::self.node_dof]
        Thetaz = Uloc[2::self.node_dof]
        Ut_def = (1-s)/2*Ut[0] + (1+s)/2*Ut[1]
        Un_def = 1/4.*(1-s)**2*(2+s)*Un[0] + L/8.*(1-s)**2*(1+s)*Thetaz[0] + \
                 1/4.*(1+s)**2*(2-s)*Un[1] + L/8.*(1+s)**2*(s-1)*Thetaz[1]
        
        x_def = x[0]*(1-s)/2.+x[1]*(1+s)/2.+Ut_def*t[0]-Un_def*t[1]
        y_def = y[0]*(1-s)/2.+y[1]*(1+s)/2.+Ut_def*t[1]+Un_def*t[0]
        return x_def,y_def
        
    def stresses(self,Ue,mat,sect):
        """ Compute generalized stresses
            
            .. math:: \{\Sigma\} = \\langle N,M_1,M_2,V\\rangle^T
        
        where :math:`M_1` is the bending moment evaluated at the first node and
        :math:`M_2` the bending moment evaluated at the second node
        
        Parameters
        ----------
        Ue : ndarray
            nodal values of the displacement
        """
        L = self.measure()
        E = mat.Young_modulus
        S = sect.area
        I = sect.inertia
        
        R = self.rotation_matrix()
        
        Uloc = np.dot(R,Ue)
        B1 = self.compute_Be_matrix(0)
        B2 = self.compute_Be_matrix(L)
        strain = np.zeros((4,))
        # axial strain and curvature at x=0
        strain[:2] = np.dot(B1,Uloc)
        # curvature at x=L
        strain[2] = np.dot(B2,Uloc)[1]
        # shear strain
        strain[3] = np.dot(np.array([0,-12./L**3,-6./L**2,0,12./L**3,-6./L**2]),Uloc)
        sig= E*np.dot(np.diag((S,I,I,I)),strain)
        return sig
