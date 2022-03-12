#coding:utf-8

import math as m
import numpy as np
import matplotlib.pyplot as plt


# Coding a part of waves and vibration course
# Exo 2 - 8

# Exo parameters
f0 = 1    # [Hz]
L  = 25*(10**-2)
k  = (f0**2)*15*(10**3)*((2*np.pi)**2)    # [kg.s-2] Stifness
w0 = (2*np.pi) * f0   # natural pulsation [rad.s-1]
Zeta = 10**-2

A = 10*(10**-2)    # Amplitude of the bump function

N = 5000   # length of frequencies
f = np.linspace(0, 2, N)
Z = A/(((2*(np.pi))**2) * pow((f0**2 - f**2)**2 + (2*Zeta*f0*f)**2,1/2))

F0 = np.zeros(N)

F0[0:N] = f0

plt.plot(f, Z, 'r', F0, Z, 'b--')

plt.xlabel('f [Hz]')
plt.ylabel('Amplitude Z [m]')
plt.title('Amplitude Z as a function of the frequency f')

VA = np.zeros(N) # VA stands for vibratory amplitude
LF = np.zeros(N)

for i in range(0,N):
    if Z[i] < (10**-2):
        VA[i] = Z[i]
        LF[i] = f[i]    # LF stands for limit frequency
        
        
for i in range(0,N):    # Sorting out the values of B as from lowest to biggest
    for j in range(0,N):
        if VA[i] < VA[j]: # We work with the frequency associated to the maximal amplitude
            c = VA[i]
            VA[i] = VA[j]
            VA[j] = c
            
            # We sort out the frequencies as well as we want to associate the Max VA to the LF
            c = LF[i]
            LF[i] = LF[j]
            LF[j] = c
        
Z_lim = VA[N-1]       
F_lim = LF[N-1]

# We then deduce the limit speed to fulfill the condition Z < 1 cm

V = L*F_lim*3.6     # Convertion to km.h-1
print(V)
print(Z_lim)

# -------------------------------------------------------

Beta = f/f0
plt.plot(Beta, np.arctan(2*Zeta*Beta/(Beta**2 - 1)), 'g', linewidth = 2)

plt.title(r'Phase \psy')
plt.xlabel('Beta [0]')
plt.ylabel('Phase [°]')

# asymptotes expected in pi/2 and -pi/2