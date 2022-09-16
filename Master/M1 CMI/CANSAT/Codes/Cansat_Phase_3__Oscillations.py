#coding:utf-8

import math as m
import numpy as np
import matplotlib.pyplot as plt


# Coding amplitude - Phase 3

# Exo arbitrary parameters
k    = 1.6 * 10**4      # [N/m]
c    = 2.5                # viscous stiffness

N    = 1000             # number of iteration (higher, the more precise)

# Fixed parameters
# force
F_0 = 9.81 - 5    # [N] Intensity oscillations

# mass
m_ressort = 5*10**-3 # [kg]
m_disque  = 10**-1   # [kg]
mass = 3*m_ressort + m_disque

# frequencies
Zeta = c/(2*m.sqrt(mass*k))

w0 = m.sqrt(k/mass)         # natural pulsation [rad.s-1]

f = np.linspace(0,200,N)
w = 2*np.pi*f

# Amplitude
X = (F_0/mass)/(pow((2*(w0**2) - w**2)**2 + (4*Zeta*w0*w)**2,1/2))

# Vector to print the frequencies f_0 & f_d
F0 = np.zeros(N)
F  = np.zeros(N)

F0[0:N] = w0/(2*np.pi)
F[0:N]  = 2*w0/(2*np.pi)*pow(1/2 - Zeta**2,1/2)

# Display results
plt.plot(f, X, 'r', label = "|X|")
plt.plot(F0, X, 'b--', label = "f_0")
plt.plot(F, X, 'g', label = "f_d")

plt.xlabel('f [Hz]')
plt.ylabel('|X(f)| [m]')
plt.title('Amplitude X as a function of the frequency f')
plt.legend()

# Save result as a graphic
# plt.savefig("3__c_3", dpi=3000)


