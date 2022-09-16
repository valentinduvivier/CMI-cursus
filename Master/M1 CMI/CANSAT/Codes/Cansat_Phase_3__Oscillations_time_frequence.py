#coding:utf-8

import math as m
import numpy as np
import matplotlib.pyplot as plt


# Coding position as a function of f and t - Phase 3

# Exo parameters
k    = 1.6 * 10**3      # [N/m]
c    = 3                # viscous stiffness

N    = 10**3            # number of iteration (higher, the more precise)

# mass
m_ressort = 5*10**-3 # [kg]
m_disque = 10*10**-3 # [kg]
mass = 3*m_ressort + m_disque

# Fixed parameters
# frequencies
Zeta = c/(2*m.sqrt(mass*k))

w0 = m.sqrt(k/mass)         # natural pulsation [rad.s-1]
wd = 2*w0*pow(1/2 - Zeta**2,1/2)

F_0 = 9.81 - 5    # [N] Intensity oscillations

# Initial parameter
x0 = 5*10**-2   # [m]     Initial condition on the position
v0 = 0          # [m.s-1] Initial condition on the speed

# mass
m_ressort = 5*10**-3 # [kg]
m_disque = 10*10**-3 # [kg]
mass = 3*m_ressort + m_disque

# Vectors
t = np.linspace(0,100,N)
f = np.linspace(0,100,N)

# Simplifications
w = 2*np.pi*f
Beta = w/w0

phi_p = np.zeros(N)
phi_d = np.zeros(N)
Xp = np.zeros(N)
X0 = np.zeros(N)


for i in range (N):
    phi_p[i] = m.atan(4*Zeta*Beta[i]/(2 - Beta[i]**2))
    Xp[i]    = (F_0/mass)/pow((2*w0**2 - w[i])**2 + (4*w0*w[i]*Zeta)**2, 1/2)
    phi_d[i] = m.atan((-2*Zeta*w0 + Xp[i]*w[i])/(Xp[i]*m.cos(phi_p[i]) - x0))
    X0[i]    = (x0 - Xp[i]*m.cos(phi_p[i]))/m.cos(phi_d[i])

x = np.zeros((N,N))

for i in range(N):
    for j in range(N):
        x[i,j] =  X0[j] * m.exp(-2*w0*t[i]) * (m.cos(wd*t[i] - phi_d[j])) + Xp[j]*m.cos(w[j]*t[i] - phi_p[j])

#for i in range(N-1):
plt.plot(t, x[1,:], label = "f = {} ".format(f[i]))

plt.xlabel('t [s]')
plt.ylabel('x [m]')

plt.title('Variation spring position over time')    # We assume the variations as being for the point x_c (center of the spring under homogeneous and isotrope conditions)
# plt.legend()

