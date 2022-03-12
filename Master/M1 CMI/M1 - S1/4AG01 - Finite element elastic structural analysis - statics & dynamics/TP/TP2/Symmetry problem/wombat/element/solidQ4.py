# -*- coding: utf-8 -*-
"""
Module for class `SolidQ4`

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
from .generic_element import *
from .trace_elements import TraceSolidT3

class SolidQ4(Quadrangle):
    """ A 2D linear quadrangular element for continuum mechanics

    """
    trace = TraceSolidT3        # defines the corresponding trace element

    def __init__(self,node_list=[],tag=1):
        """
        Parameters
        ----------

        node_list : list
            list containing two nodes
        tag : int,str
            tag of physical group
        """
        Quadrangle.__init__(self,node_list,tag)
        self.elem_type = 'Q4'
        self.el_dof = 8
        self.node_dof = 2
        self.ngauss = 4
        self.nb_stresses = 4*self.ngauss  # 4 stresses at ngauss Gauss points in element

        self.ag,self.wg = self.gauss_quadrature(self.ngauss)
        self.detJ = [self.jacobian(a[0],a[1])[0] for a in self.ag]
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
        N : ndarray shape = ?
            array of shape functions :math:`[N]` evaluated at :math:`(\\xi,\\eta)`
        DN : ndarray shape = ?
            array of shape functions derivatives :math:`[\\nabla N]` evaluated at :math:`(\\xi,\\eta)`
        """
        N = np.array([(1-xi)*(1-eta)/4.,
                      (1+xi)*(1-eta)/4.,
                      (1+xi)*(1+eta)/4.,
                      (1-xi)*(1+eta)/4.])
        DN = np.array([[-(1-eta)/4.,
                        (1-eta)/4.,
                        (1+eta)/4.,
                        -(1+eta)/4.],[
                        -(1-xi)/4.,
                        -(1+xi)/4.,
                        (1+xi)/4.,
                        (1-xi)/4.]])
        return N,DN

    def jacobian(self,xi,eta):
        """ Computes quantities related to the jacobian of the element

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
        N, DN = self.shape_functions(xi,eta)
        Jac = np.dot(DN,T)
        detJ = np.linalg.det(Jac)
        assert detJ>0, "Jacobian of element %i is negative." % self._id
        invJac = np.linalg.inv(Jac)
        return detJ,Jac,invJac

    def compute_Be_matrix(self,xi,eta):
        """ Local strain matrix :math:`[B_e]` such that

        .. math:: [B_e]\{U_e\} = \\begin{Bmatrix} \\varepsilon_{xx} \\\\ \\varepsilon_{yy} \\\\ 2\\varepsilon_{xy} \end{Bmatrix}

        evaluated at :math:`(\\xi,\\eta)`, shape = (3,12)
        """
        N,DN = self.shape_functions(xi,eta)
        detJ, Jac, invJac = self.jacobian(xi,eta)
        GN = np.dot(invJac,DN)
        Be = np.zeros((3,8))
        Be[0,0::2] = GN[0,:]
        Be[1,1::2] = GN[1,:]
        Be[2,0::2] = GN[1,:]
        Be[2,1::2] = GN[0,:]
        return Be


    def elementary_stiffness(self,mat,sect=None):
        """ Elementary stiffness matrix :math:`[K_e]` shape=(12,12)"""
        C = mat.C_matrix()
        C_plan = np.copy(C)
        C_plan[2,2] = 0
        C_cis = np.zeros((3,3))
        C_cis[2,2] = C[2,2]
        Ke = np.zeros((8,8))
        for i in range(self.ngauss):
#            Ke += self.detJ[i]*self.wg[i]*np.dot(self.Be[i].T,np.dot(C,self.Be[i]))
            Ke += self.detJ[i]*self.wg[i]*np.dot(self.Be[i].T,np.dot(C_plan,self.Be[i]))
        detJ0 = self.jacobian(0.,0.)[0]
        Be0 = self.compute_Be_matrix(0.,0.)
        Ke += 4*detJ0*np.dot(Be0.T,np.dot(C_cis,Be0))
        return Ke


    def elementary_mass(self,mat,sect=None,lumped=False):
        pass

    def elementary_thermal_vector(self,mat,dilat):
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
        fe = np.zeros((8,))
        for i in range(self.ngauss):
            N,DN = self.shape_functions(self.ag[i,0],self.ag[i,1])
            fe += self.detJ[i]*self.wg[i]*np.kron(N,np.array([fx,fy]))
        return fe

    def deformation(self,Ue):
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

        x = self.node_coor()[:,0]
        y = self.node_coor()[:,1]
        x_def = x + Ux
        y_def = y + Uy
        x_def = np.append(x_def, x_def[0])
        y_def = np.append(y_def, y_def[0])
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
        C = mat.C_matrix()
        sig_g = np.zeros((4*self.ngauss,))
        for i in range(self.ngauss):
            Eps = np.dot(self.Be[i],Ue)
            sig_g[[4*i,4*i+1,4*i+3]] =  np.dot(C,Eps)
            sig_g[4*i+2] = mat.compute_sigzz(Eps)
        return sig_g


