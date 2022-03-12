# -*- coding: utf-8 -*-
"""
Module for class `SolidT6`

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
from .generic_element import *
from .trace_elements import TraceSolidT6

class SolidT6(Triangle6):
    """ A 2D quadratic triangular element for continuum mechanics
    
        SolidT6 is a Triangle6-type element (6 nodes) in 2D
        with 2 degrees of freedom/node and with **straight sides** (affine interpolation of the geometry):\n 
        
        - **Kinematics**:  horizontal, vertical displacement \
             with a quadratic interpolation inside the element 
        
        .. math:: \{U\}=\\langle u_x^1,u_y^1,\\ldots,u_x^6,u_y^6\\rangle^T 
            
        - **Strains** :  plane components of :math:`\\underline{\\underline{\\varepsilon}} \
           = (\\underline{\\nabla u} + \\underline{\\nabla u}^T)/2`
           :math:`\\varepsilon_{xx},\\varepsilon_{yy},\\varepsilon_{xy}` (linear)
        
        - **Stresses**: :math:`\\{\\sigma_{xx},\\sigma_{yy},\\sigma_{zz},\\sigma_{xy}\\}` (linear)
        
            .. note:: the out-of-plane stress :math:`\\sigma_{zz}` is not necesseary for building the stiffness \
            matrix but is still given as an output (computed from :math:`\\sigma_{zz}= \
            \\lambda(\\varepsilon_{xx}+\\varepsilon_{yy})`)
    """
    trace = TraceSolidT6        # defines the corresponding trace element
    
    def __init__(self,node_list=[],tag=1):  
        """
        Parameters
        ----------
        
        node_list : list
            list containing two nodes
        tag : int,str
            tag of physical group
        """      
        Triangle6.__init__(self,node_list,tag)
        self.elem_type = 'T6'
        self.el_dof = 12
        self.node_dof = 2
        self.ngauss = 3
        self.nb_stresses = 4*self.ngauss  # 4 stresses at ngauss Gauss points in element
        
        self.kin_field_names = ['U_x','U_y']
        self.strain_field_names = ['eps_xx','eps_yy','2eps_xy']
        self.int_forces_field_names = ['sig_xx','sig_yy','sig_zz','sig_xy']
        self.ext_forces_field_names = ['F_x','F_y']
        
        self.detJ,self.Jac,self.invJac = self.jacobian()
        self.A = self.detJ/2.
        self.ag,self.wg = self.gauss_quadrature(self.ngauss)
        self.Be = [self.compute_Be_matrix(a[0],a[1]) for a in self.ag]
            
    def shape_functions(self,xi,eta):
        """ Returns the shape functions and its derivatives
        
        Parameters
        -----------
        xi : float
            coordinate of point :math:`\\xi` in the reference space, belongs to :math:`[0;1]`   
        eta : float
            coordinate of point :math:`\\eta` in the reference space, belongs to :math:`[0;1]`
        
        Returns
        --------
        N : ndarray shape = (6,)
            array of shape functions :math:`[N]` evaluated at :math:`(\\xi,\\eta)`
        DN : ndarray shape = (2,6)
            array of shape functions derivatives :math:`[\\nabla N]` evaluated at :math:`(\\xi,\\eta)`
        """
        #TO BE CODED

    def jacobian(self):
        """ Computes quantities related to the jacobian of the element (here constant)

        Returns
        --------
        detJ : float
            determinant of the jacobian matrix (must be strictly positive)
        Jac : ndarray shape (2,2)
            jacobian matrix :math:`[J]`
        invJac : ndarray shape (2,2)
            inverse of the jacobian matrix :math:`[J]^{-1}`
        """
        T = self.node_coor()
        Jac = np.array([T[0,:]-T[2,:],T[1,:]-T[2,:]])
        detJ = np.linalg.det(Jac)
        assert detJ>0, "Jacobian of element %i is negative." % self._id
        invJac = np.linalg.inv(Jac)
        return detJ,Jac,invJac
        
    def compute_Be_matrix(self,xi,eta):
        """ Local strain matrix :math:`[B_e]` such that
        
        .. math:: [B_e]\{U_e\} = \\begin{Bmatrix} \\varepsilon_{xx} \\\\ \\varepsilon_{yy} \\\\ 2\\varepsilon_{xy} \end{Bmatrix} 
        
        evaluated at :math:`(\\xi,\\eta)`, shape = (3,12)
        """
        #TO BE CODED
            
    
    def elementary_stiffness(self,mat,sect=None):     
        """ Elementary stiffness matrix :math:`[K_e]` shape=(12,12)"""    
        #TO BE CODED
        
    def elementary_thermal_vector(self,mat,dilat):
        # TODO
        pass
        
    
    def elementary_distributed_forces(self,el_force):
        """ Elementary force vector for uniform distributed loading
        
        Parameters
        ----------
        el_force = [fx,fy,cz] : array,list
            contains uniformly distributed forces :math:`(f_x,f_y)`
            
           .. note:: for the SolidT6 element cz is ignored           
        """
        fx,fy,cz = el_force
        fe = np.zeros((12,))
        for i in range(self.ngauss):
            N,DN = self.shape_functions(self.ag[i,0],self.ag[i,1])
            fe += self.detJ*self.wg[i]*np.kron(N,np.array([fx,fy]))
        return fe
    
    
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
            returns deformed position of m points along the element boundary
        """
        Ux = Ue[0::self.node_dof]
        Uy = Ue[1::self.node_dof]
        s = np.linspace(-1,1,m)
        # repeat twice along axis 0
        s = np.tile(s,(1,1)).T
        Nq1 = s*(s-1)/2
        Nq2 = s*(1+s)/2
        Nq3 = 1-s**2
        z = 0*s
        S1 = np.concatenate((Nq1,Nq2,z,Nq3,z,z),axis=1)
        S2 = np.concatenate((z,Nq1,Nq2,z,Nq3,z),axis=1)
        S3 = np.concatenate((Nq2,z,Nq1,z,z,Nq3),axis=1)
        S = np.concatenate((S1,S2,S3),axis=0)
        x = self.node_coor()[:,0]
        y = self.node_coor()[:,1]
        x_def = np.dot(S,x+Ux)
        y_def = np.dot(S,y+Uy)
        return x_def,y_def

    def stresses(self,Ue,mat,sect=None):
        """ Compute stress state evaluated at Gauss points, shape = (4*ngauss,)
            
            .. math:: \{\Sigma\} = \\begin{Bmatrix} \\Sigma^1 \\\\ \\vdots \\\\ \\Sigma^{ngauss} \end{Bmatrix} \\text{ with } \{\Sigma^g\} = \\begin{Bmatrix} \\sigma_{xx}^g \\\\ \\sigma_{yy}^g \\\\ \\sigma_{zz}^g \\\\ \\sigma_{xy}^g \\end{Bmatrix}

        .. note:: :math:`\\sigma_{zz}` is not used for the computation but computed from the strain       
        
        Parameters
        ----------
        Ue : ndarray
            nodal values of the displacement
        """      
        #TO BE CODED

    def internal_forces(self,Sige):
        """ Returns elemental contribution of a stress state Sige to internal forces vector Fint """
        fe = np.zeros((12,))
        for i in range(self.ngauss):
            fe += self.wg[i]*self.detJ*np.dot(self.Be[i].T,Sige[[0,1,3]])
        return fe

