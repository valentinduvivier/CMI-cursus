#coding:utf-8

import math as m
import numpy as np
import matplotlib.pyplot as plt


# Coding angular position - Phase 4

"""
    I measured from Catia model
"""

# Exo arbitrary parameters
K    = 30       # [N.m.rad-1]
C    = 3*10**-4   # viscous stiffness
C_1  = 3*10**-3

theta_ini = np.pi/6       # Angle between the vertical and the rod once landing
I         = 4.04*10**-5

N         = 200     # number of iteration (higher, the more precise)

# Fixed parameters
# Dimensions
v_choc      = 7             # Speed of the module when landing
L           = 10*10**-2     # Position center of gravity of a rod

# Forces - Momentum
F_reaction  = 9.81 - 5            # Force applied by the floor on the system as a reaction
M           = L * F_reaction * m.sin(theta_ini)   # Moment due to the reaction force

# Frequency
wd = pow(4*K**2 - C**2, 1/2)/(2*I)
wd_1 = pow(4*K**2 - C_1**2, 1/2)/(2*I)

phi_d  = m.atan(-C/(2*I*wd))
phi_d_1  = m.atan(-C/(2*I*wd))

Theta0 = (theta_ini - M/K)/m.cos(phi_d) 
Theta0_1 = (theta_ini - M/K)/m.cos(phi_d_1) 

# Vector time
t = np.linspace(0,1,N)

# Position
theta  = np.zeros(N)
theta_1= np.zeros(N)
dtheta = np.zeros(N)

for i in range(N):
    # Position
    theta[i] =  Theta0 * m.exp(-C*t[i]/(2*I)) * m.cos(wd*t[i] - phi_d) + M/K
    theta_1[i] =  Theta0 * m.exp(-C_1*t[i]/(2*I)) * m.cos(wd_1*t[i] - phi_d_1) + M/K
    # Speed
    # dtheta[i] =  ( - C/(2*I) * Theta0 * m.exp(-C*t[i]/(2*I)) * m.cos(wd*t[i] - phi_d) ) + ( - Theta0 * wd * m.exp(-C*t[i]/(2*I)) * m.sin(wd*t[i] - phi_d) )

# Displaying results
# Position
plt.plot(t, theta*360/(2*np.pi), 'g', label = "$\Theta_{ini}$ = $\pi$/6 & C << 1")
plt.plot(t, theta_1*360/(2*np.pi), 'r', label = "$\Theta_{ini}$ = $\pi$/6 & C >> 1")

# Speed
# plt.plot(t, dtheta, 'g', label = "Zeta < 0.707")


plt.xlabel('t [s]')
plt.ylabel('theta [°]')
plt.title('Variation spring angular position over time')    # We assume the variations as being for the point x_c (center of the spring under homogeneous and isotrope conditions)
plt.legend()

# Saving results
# plt.savefig("4__c_001.png", dpi=3000)


