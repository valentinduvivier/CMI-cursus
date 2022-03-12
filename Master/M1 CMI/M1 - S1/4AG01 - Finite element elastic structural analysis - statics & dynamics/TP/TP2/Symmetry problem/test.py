# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 15:25:19 2020

@author: vltn0
"""

import numpy as np
import matplotlib.pyplot as plt

Length = 40*10**-2      # [m]

A = np.array([
             [1, 0.568],
             [2, 0.496],
             [3, 0.462],
             [4, 1.497],
             [5, 2.009],
             [6, 2.724],
             [7, 3.329],
             [8, 35.03]
             ])

sigma = np.array([
                 [],
                 [],
                 [],
                 [],
                 [],
                 [],
                 [],
                 [],
                 []                 
                 ])

A = A*Length/2

plt.plot(A[:,0], A[:,1], 'r-')

plt.title(r'Stress concentration factor - $K_G$')
plt.xlabel(r'n --- $a=n*b$')
plt.ylabel('$K_G$')

plt.grid('True')