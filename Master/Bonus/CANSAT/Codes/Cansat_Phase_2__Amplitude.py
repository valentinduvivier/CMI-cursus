#coding:utf-8

import math as m
import numpy as np
import matplotlib.pyplot as plt


# Coding amplitude - Phase 2

# Exo arbitrary parameters
k    = 1.6 * 10**3       # [N/m]
c    = 3        # viscous stiffness

N = 1000        # number of iteration (higher, the more precise)

# Fixed parameters
# force
F_choc = 9.81    # [N] Force absorbed during the choc

# mass
m_ressort = 5*10**-3    # [kg]
m_disque  = 10**-1   # [kg]
mass      = 3*m_ressort + m_disque

# position
x0        = 5*10**-2    # [m]

# Frequencies
Zeta = c/(2*m.sqrt(mass*k))
w0 = m.sqrt(k/mass)  # natural pulsation [rad.s-1]

f = np.linspace(0,400,N)   # Force's associated frequency
w = 2*np.pi*f

# Amplitude
X = pow((((F_choc/mass + 4*Zeta*w0*x0)**2 + (w*x0)**2))/((2*(w0**2) - w**2)**2 + (4*Zeta*w0*w)**2), 1/2)

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

# Save results
# plt.savefig("2__c_3.png", dpi = 3000)

