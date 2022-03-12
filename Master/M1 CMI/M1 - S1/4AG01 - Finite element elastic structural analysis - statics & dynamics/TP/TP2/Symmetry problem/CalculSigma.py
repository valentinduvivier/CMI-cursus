# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 17:11:03 2020

@author: vltn0
"""

import numpy as np

def sig_rr(f, a, r, theta):
    ratio = a/r
    return f/2*((1 - ratio**2) + np.cos(2*(theta - np.pi/2))*(1 - 4*ratio**2 + 3*ratio**4))

def sig_tt(f, a, r, theta):
    ratio = a/r
    return f/2*((1 + ratio**2) - np.cos(2*(theta - np.pi/2))*(1 + 3*ratio**4))

def sig_rt(f, a, r, theta):
    ratio = a/r
    return -f/2*(1 + 2*ratio**2 - 3*ratio**4)*np.sin(2*(theta - np.pi/2))
