#coding:utf-8

import math as m
import numpy as np
import matplotlib.pyplot as plt


# Coding position - Phase 1

g = 9.81

# Exo arbitrary parameters
k    = 1.6 * 10**3       # [N/m]
c    = 3        # viscous stiffness

N = 1000        # number of iteration (higher, the more precise)

# Fixed parameters
x0 = 5*10**-2   # Initial condition on the position
v0 = 0          # Initial condition on the speed

m_ressort = 5*10**-3 # [kg]
m_disque  = 10**-1   # [kg]
mass = 3*m_ressort + m_disque

# Frequencies
Zeta = c/(2*m.sqrt(mass*k))

w0 = m.sqrt(k/mass)  # natural pulsation [rad.s-1]
wd = 2*w0*pow(1/2 - Zeta**2,1/2)

# Amplitude
X0 = pow(x0**2 + (v0/(2*w0*wd*Zeta))**2, 1/2)

# Vector time
t = np.linspace(0,0.2,N)

# Position
x = np.zeros(N)
for i in range(N):
    x[i] =  X0 * m.exp(-2*w0*t[i]) * (m.cos(wd*t[i] - m.atan(-v0/(2*w0*wd*Zeta*x0)))) + mass*g/k

# Displaying results
plt.plot(t, x, 'g', label = "système amorti - $\zeta$ < 1")

plt.xlabel('t [s]')
plt.ylabel('x [m]')
plt.title('Variation spring position over time')    # We assume the variations as being for the point x_c (center of the spring under homogeneous and isotrope conditions)
plt.legend()

# Saving results
# plt.savefig("1__c_3", dpi=3000)


