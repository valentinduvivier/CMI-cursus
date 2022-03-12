# -*- coding: utf-8 -*-
"""
Module for class `Bar2D`

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
#spatinet change
from .generic_element import *
#spatinet change

class Bar2D(Segment): 
    """ A 2D truss element
    
        Bar2D is a Segment-type element (2 end nodes) in 2D
        with 2 degrees of freedom/node :\n 
        
        - **Kinematics**: horizontal, vertical displacement \
        with a linear interpolation inside the element 
        
            .. math:: \{U\}=\\langle u_x^1,u_y^1,u_x^2,u_y^2\\rangle^T
            
        - **Strains**: axial strain :math:`\epsilon`  (constant)
        - **Stresses**: normal force :math:`N` (constant)
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
        self.el_dof = 4
        self.node_dof = 2
        self.nb_stresses = 1
        self.kin_field_names = ['U_x','U_y']
        self.strain_field_names = ['eps']
        self.stresses_field_names = ['N']
        self.ext_forces_field_names = ['F_x','F_y']
    
    def rotation_matrix(self):
        """
        Rotation matrix :math:`[R]` from global to local frame 
        
        shape = (2,4)
        """
        T = self.nodes.coor
        tang = T[1,:]-T[0,:]
        L = self.measure()
        t = tang/L
        
        # rotation matrix [r]{U_x,U_y} = {U_t}
        r = np.array([t[0],t[1]])
        # R is such that [R]{U_x1,U_y1,U_x2,U_y2} = {U_t1,U_t2}
        return np.kron(np.eye(2),r) # get a block matrix by repeating [r] twice along the diagonal, 
                                 # equivalent to R = np.bmat([[r,np.zeros((1,2))],[np.zeros((1,2)),r]])
        
    def elementary_stiffness(self,mat,sect):
        """ Elementary stiffness matrix :math:`[K_e]` shape=(4,4)
        
        elementary stiffness in local frame is
        
        .. math:: [K_{e,loc}]=\\dfrac{ES}{L}\\begin{bmatrix} 1 & -1 \\\\ -1 & 1\end{bmatrix}
        """ 
        L = self.measure()
        E = mat.Young_modulus
        S = sect.area
        R = self.rotation_matrix()
        
        # elementary stiffness matrix in local frame
        Ke_loc = E*S/L*(np.array([[1,-1],[-1,1]]))
        
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
        
        return E*S*dilat*np.hstack((-tang,tang))
        
    def elementary_distributed_forces(self,el_force):
        """ Elementary force vector for uniform distributed loading
        
        Parameters
        ----------
        el_force = [fx,fy,cz] : array,list
            contains uniformly distributed forces :math:`(f_x,f_y)`
            
           .. note:: for the Bar2D element cz is ignored
        """
        L = self.measure()
        fx,fy,cz = el_force
        return np.tile([fx,fy],2)*L/2.
    
    def deformation(self,Ue):
        """ Interpolation of the deformed element
        
        Parameters
        ----------
        Ue : ndarray
            nodal displacement of the current element
        
        Returns
        -------
        x_def,y_def : ndarray
            returns deformed position of the two end nodes
        """
        Ux = Ue[0::self.node_dof]
        Uy = Ue[1::self.node_dof]
        s = np.linspace(-1,1,2)
        x = self.node_coor()[:,0]
        y = self.node_coor()[:,1]
        x_def = (1-s)/2.*(Ux[0]+x[0])+(1+s)/2.*(Ux[1]+x[1])
        y_def = (1-s)/2.*(Uy[0]+y[0])+(1+s)/2.*(Uy[1]+y[1])
        return x_def,y_def
        
    def stresses(self,Ue,mat,sect):
        """ Compute generalized stresses
            
            .. math:: \{\Sigma\} = \{N\}
        
        Parameters
        ----------
        Ue : ndarray
            nodal values of the displacement
        """
        L = self.measure()
        E = mat.Young_modulus
        S = sect.area
        
        #tangential displacement
        Ut = np.dot(self.rotation_matrix(),Ue)
        
        return E*S/L*(Ut[1]-Ut[0])
        