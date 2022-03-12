# -*- coding: utf-8 -*-
"""
Module for generating and reading Gmsh mesh files

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
from .node import Node2D, NodeGroup
from .element import ElementGroup, Segment,Segment3,Triangle, Triangle6
from .mesh import Mesh
import numpy as np
import os
import meshio
from sys import platform

def call_gmsh(geofile,elem_type):
    """ Calls ``Gmsh`` to generate a mesh
    
    Parameters
    ----------
    geofile : str
        path to GEO file, must end in ".geo"
    elem_type : {:class:`SolidT3`,:class:`SolidT6`}
        type of generated elements
    
    Returns
    -------
    regions
        a :class:`Mesh <mesh.Mesh>` object 
    """
    if issubclass(elem_type,Triangle6):
        order = 2
    else:
        order = 1

    # *** HERE ADD YOUR PATH TO GMSH
    #command = 'your/gmsh/folder/gmsh -order '+str(order)+' '+geofile+' -2'  
    
    # Example (RC)
    command = '/Applications/Gmsh.app/Contents/MacOs/gmsh -order '+str(order)+' '+geofile+' -2'  
  
    os.system(command)
    print ("Finished meshing!\n")
    mshfile = geofile[:-3]+'msh'
    return read_msh(mshfile,elem_type)
    
def read_msh(mshfile,elem_type):
    """Reads ``Gmsh`` 2.0 mesh files using ``meshio`` package
    
    Parameters
    ----------
    mshfile : str
        path to MSH file, must end in ".msh"
    
    Returns
    -------
    regions
        a :class:`Mesh <mesh.Mesh>` like object with abstract geometrical entities (Segment, Triangles, etc.) containing
        mesh cells, facets and nodes
    """

    # R. Cornaggia : the version below corresponds to the mesh data structure
    # given by meshio v. 4.3.5
    # Ask for help if you have a former version that does not works
    
    mesh = meshio.read(mshfile)    
    cells = mesh.cells_dict
    cell_data = mesh.cell_data_dict["gmsh:physical"]
    field_data = mesh.field_data    
    
    region_names = {}
    for keys, values in field_data.items():
        region_names.update({values[0]: keys})
    def get_name(v):
        name = region_names[v]
        if name is not None:
            return name
        else:
            return v      

    Nno = mesh.points.shape[0]

    nodes = NodeGroup([Node2D(mesh.points[i,:2]) for i in range(Nno)])
    if "line" in cells:
        cell_l = cells["line"]
        cell_data_l = cell_data["line"]
    else:
        cell_l = np.empty(shape=(0,0))
    if "line3" in cells:
        cell_l3 = cells["line3"]
        cell_data_l3 = cell_data["line3"]
    else:
        cell_l3 = np.empty(shape=(0,0))
    if "triangle" in cells:
        cell_T3 = cells["triangle"]
        cell_data_T3 = cell_data["triangle"]
    else:
        cell_T3 = np.empty(shape=(0,0))
    if"triangle6" in cells:
        cell_T6 = cells["triangle6"]
        cell_data_T6 = cell_data["triangle6"]
    else:
        cell_T6 = np.empty(shape=(0,0))

# ************************************ end changes COrnaggia

        
    if issubclass(elem_type,Triangle) and not issubclass(elem_type,Triangle6):
        if (len(cell_T6)>0 or len(cell_l3)>0):
            raise(ValueError,"Trying to use 3-noded elements with 6-noded triangles or 3-noded lines in .msh file")
        else:
            elem_list = [elem_type([nodes.node_list[cell_T3[i,j]] for j in range(3)], tag=get_name(cell_data_T3[i])) for i in range(cell_T3.shape[0])]
            bd_list = [elem_type.trace([nodes.node_list[cell_l[i,j]] for j in range(2)], tag=get_name(cell_data_l[i])) for i in range(cell_l.shape[0])]                      
    elif issubclass(elem_type,Triangle6):
        if len(cell_T3)>0 or len(cell_l)>0:
            raise(ValueError,"Trying to use 6-nodes elements with 3-noded triangles or 2-noded lines in .msh file")
        else:
            elem_list = [elem_type([nodes.node_list[cell_T6[i,j]] for j in range(6)], tag=get_name(cell_data_T6[i])) for i in range(cell_T6.shape[0])]
            bd_list = [elem_type.trace([nodes.node_list[cell_l3[i,j]] for j in range(3)], tag=get_name(cell_data_l3[i])) for i in range(cell_l3.shape[0])]
    mesh = Mesh(elem_list)
    boundary = ElementGroup(bd_list)

    return mesh,boundary
    
