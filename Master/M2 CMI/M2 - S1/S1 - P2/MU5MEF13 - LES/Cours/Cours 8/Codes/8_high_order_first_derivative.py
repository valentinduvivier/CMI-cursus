# Script to plot derivative modified wavenumber of the form (compact/explicit):
#
# beta*f'_{i-2} + alpha*f'_{i-1} + f'_i + alpha*f'_{i+1} + beta*'_{i+2}  =
# + a*( f_{i+1} - f_{i-1} ) + b*( f_{i+2} - f_{i-2} ) + c*( f_{i+3} - f_{i-3} )
# + d*( f_{i+4} - f_{i-4} ) + e*( f_{i+5} - f_{i-5} ) + ..
#
# This works only for centered schemes.
import numpy as np
import matplotlib.pyplot as plt

npt = 200                           # number of points
kdx = np.linspace( 0., np.pi, npt ) # reduced wavenumber

# Function for modified wavenumber computation
def modified_wavenumber( scheme ):
   "Derivative transfer function"
   num = 0
   for k in range(len(scheme['exp_coeff'])):
      num = num + 2.*scheme['exp_coeff'][k]*np.sin( (k+1)*kdx )
   den = 1
   for k in range(len(scheme['imp_coeff'])):
      den = den + 2.*scheme['imp_coeff'][k]*np.cos( (k+1)*kdx )
   return num / den

list_schemes = list() # List of FD schemes to plot

# ==============================================================================
#                       Explicit Standard Centred Differences
# ==============================================================================
# coefficients 2nd order explicit 3 pts
list_schemes.append( { 'exp_coeff': [1/2]
                     , 'imp_coeff': []
                     , 'style'    : ['magenta', 'solid', '', 'STD 3pt o2']})
# coefficients 4th order explicit 5 pts
list_schemes.append( { 'exp_coeff': [2/3, -1/12]
                     , 'imp_coeff': []
                     , 'style'    : ['green', 'solid', '', 'STD 5pt o4']})
# coefficients 6th order explicit 7 pts
list_schemes.append( { 'exp_coeff': [3/4, -3/20, 1/60]
                     , 'imp_coeff': []
                     , 'style'    : ['orange', 'solid', '', 'STD 7pt o6']})
# coefficients 8th order explicit 9 pts
list_schemes.append( { 'exp_coeff': [4/5, -1/5, 4/105, -1/280]
                     , 'imp_coeff': []
                     , 'style'    : ['blue', 'solid', '', 'STD 9pt o8']})
# coefficients 10th order explicit 11 pts
list_schemes.append( { 'exp_coeff': [5/6, -5/21, 5/84, -5/504, 1/1260]
                     , 'imp_coeff': []
                     , 'style'    : ['red', 'solid', '', 'STD 11pt o10']})
# coefficients 12th order explicit 13 pts
list_schemes.append( { 'exp_coeff': [6/7, -15/56, 5/63, -1/56, 3/1155, -1/5544]
                     , 'imp_coeff': []
                     , 'style'    : ['gray', 'solid', '', 'STD 13pt o12']})
# ==============================================================================
#                  Implicit Centred Differences (Lele JPC 1992)
# ==============================================================================
# 4th order family compact 3imp-3exp, classical pade' 4th order
alpha=1/4; beta=0
list_schemes.append( { 'exp_coeff': [1/3*(alpha+2), 1/12*(4*alpha-1)]
                     , 'imp_coeff': [alpha, beta]
                     , 'style'    : ['orange', 'dashed', '', 'CMP 3-3pt o4']})
# 4th order family compact 3imp-3exp, 6th order minimal stencil
alpha=1/3; beta=0
list_schemes.append( { 'exp_coeff': [1/3*(alpha+2), 1/12*(4*alpha-1)]
                     , 'imp_coeff': [alpha, beta]
                     , 'style'    : ['blue', 'dashed', '', 'CMP 3-3pt o6']})
# 6th order family compact 3imp-5exp, 8th order minimal stencil
alpha=3/8; beta=0
list_schemes.append( { 'exp_coeff': [1/12*(alpha+9), 1/60*(32*alpha-9), 1/60*(-3*alpha+1)]
                     , 'imp_coeff': [alpha, beta]
                     , 'style'    : ['red', 'dashed', '', 'CMP 3-5pt o8']})
# 8th order family compact 5imp-5exp --> 10th order minimal stencil
alpha=1/2; beta=1/20*(-3+8*alpha)
list_schemes.append( { 'exp_coeff': [1/12*(12-7*alpha), 1/600*(568*alpha-183), 1/300*(9*alpha-4)]
                     , 'imp_coeff': [alpha, beta]
                     , 'style'    : ['gray', 'dashed', '', 'CMP 5-5pt o10']})
# ==============================================================================
#                       Optimized Centred Differences
# ==============================================================================
# explicit DRP 7 pts (Bogey & Bailly 2004 JCP)
list_schemes.append( { 'exp_coeff': [0.790803914666667,-0.182643131733333, 0.0248274496]
                     , 'imp_coeff': []
                     , 'style'    : ['orange', 'dotted', '', 'DRP 7pt o4']})
# explicit4DRP 9 pts (Bogey & Bailly 2004 JCP)
list_schemes.append( { 'exp_coeff': [0.841570125482,-0.244678631765, 0.059463584768, -0.007650904064]
                     , 'imp_coeff': []
                     , 'style'    : ['blue', 'dotted', '', 'DRP 9pt o4']})
# explicit DRP 11 pts (Bogey & Bailly 2004 JCP)
list_schemes.append( { 'exp_coeff': [0.872756993962,-0.286511173973, 0.090320001280, -0.020779405824, 0.002484594688]
                     , 'imp_coeff': []
                     , 'style'    : ['red', 'dotted', '', 'DRP 11pt o4']})
# explicit DRP 13 pts (Bogey & Bailly 2004 JCP)
list_schemes.append( { 'exp_coeff': [0.907646591371,-0.337048393268, 0.133442885327, -0.045246480208, 0.011169294114, -0.001456501759]
                     , 'imp_coeff': []
                     , 'style'    : ['gray', 'dotted', '', 'DRP 13pt o4']})

#===============================================================================
#        Compute the modified wavenumber for each scheme to plot
#===============================================================================
for scheme in list_schemes:
   kstardx = modified_wavenumber(scheme)
   scheme['kstardx'] = kstardx                     # Modified wavenumber
   scheme['dispers'] = np.abs(kdx - kstardx)/np.pi # Dispersion error
   scheme['phase']   = kstardx/kdx                 # Phase error

#===============================================================================
#                             And now plot
#===============================================================================
plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')
plt.rcParams['xtick.major.pad']='8'

#===============================================================================
# Plot modified wavenumber in linear scale
# ========================================
fig, ax = plt.subplots(figsize=(8,6))
ax.plot( kdx, kdx, color='black', label='Exact' )  # Exact solution

for scheme in list_schemes:
   ax.plot(kdx, scheme['kstardx'],     color=scheme['style'][0]
                                 , linestyle=scheme['style'][1]
                                 ,    marker=scheme['style'][2]
                                 ,     label=scheme['style'][3] )
# Style formatting
ax.legend(loc=2, fontsize='x-large')
ax.set_ylim([0,np.pi])
ax.set_xlim([0,np.pi])
ax.set_xlabel(r'$k\Delta x$', fontsize=24)
ax.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
ax.set_xticklabels(['$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$',
                     r'$\frac{3\pi}{4}$', r'$\pi$'], fontsize=24)
ax.set_ylabel(r'$k^*\Delta x$', fontsize=24)
ax.set_yticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
ax.set_yticklabels([r'$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$',
                    r'$\frac{3\pi}{4}$', r'$\pi$'], fontsize=24)

plt.tight_layout()
plt.savefig('first_der_modified_wavenumber.pdf')

#===============================================================================
# Plot Dispersion error
# =====================
fig, ax = plt.subplots(figsize=(8,6))

for scheme in list_schemes:
   ax.plot(kdx, scheme['dispers'],     color=scheme['style'][0]
                                 , linestyle=scheme['style'][1]
                                 ,    marker=scheme['style'][2]
                                 ,     label=scheme['style'][3] )
# Style formatting
ax.legend(loc=2, fontsize='x-large')
ax.set_xlim(np.pi/16, np.pi)
ax.set_ylim(1e-5, 1)
ax.loglog()
ax.set_xlabel(r'$k\Delta x$', fontsize=24)
ax.set_xticks([np.pi/16, np.pi/8, np.pi/4, np.pi/2, np.pi])
ax.set_xticklabels([r'$\frac{\pi}{16}$', r'$\frac{\pi}{8}$', r'$\frac{\pi}{4}$',
                    r'$\frac{\pi}{2}$', r'$\pi$'], fontsize=24)
ax.set_ylabel(r'$|k\Delta x - k^*\Delta x|/\pi$', fontsize=24)
ax.set_yticks([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
ax.set_yticklabels([r'$10^{-6}$', r'$10^{-5}$', r'$10^{-4}$', r'$10^{-3}$',
                    r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$',], fontsize=24)
plt.tight_layout()
plt.savefig('first_der_dispersion_error.pdf')

# #===============================================================================
# Plot phase error
# ================
fig, ax = plt.subplots(figsize=(8,6))
ax.plot( kdx, kdx/kdx, color='black', label='Exact' )  # Exact solution

for scheme in list_schemes:
   ax.plot(kdx, scheme['phase'],     color=scheme['style'][0]
                               , linestyle=scheme['style'][1]
                               ,    marker=scheme['style'][2]
                               ,     label=scheme['style'][3] )
# Style formatting
ax.legend(loc=3, fontsize='x-large')
ax.set_ylim([0,1.05])
ax.set_xlim([0,np.pi])
ax.set_xlabel(r'$k\Delta x$', fontsize=24)
ax.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
ax.set_xticklabels(['$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$',
                     r'$\frac{3\pi}{4}$', r'$\pi$'], fontsize=24)
ax.set_ylabel(r'$k^*\Delta x / k\Delta x$', fontsize=24)

plt.tight_layout()
plt.savefig('first_der_phase_error.pdf')

# Show everything
plt.show()
