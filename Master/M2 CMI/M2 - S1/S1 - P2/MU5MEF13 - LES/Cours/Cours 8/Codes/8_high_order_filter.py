# Script to plot derivative modified wavenumber of the form (compact/explicit):
#
# beta*f'_{i-2} + alpha*f'_{i-1} + f'_i + alpha*f'_{i+1} + beta*'_{i+2}  =
# + a*f_{i} + b/2*( f_{i+1} - f_{i-1} ) + c/2*( f_{i+2} - f_{i-2} )
#           + d/2*( f_{i+3} - f_{i-3} ) + e/2*( f_{i+4} - f_{i-4} ) + ..
#
# This works only for centered schemes.
import numpy as np
import matplotlib.pyplot as plt

npt = 200                           # number of points
kdx = np.linspace( 0., np.pi, npt ) # reduced wavenumber

# Function for modified wavenumber computation
def transfer_function( scheme ):
   "Filter transfer function"
   num = scheme['exp_coeff'][0]
   for k in range(len(scheme['exp_coeff'])-1):
      num = num + 2*scheme['exp_coeff'][k+1]*np.cos( (k+1)*kdx )
   den = 1
   for k in range(len(scheme['imp_coeff'])):
      den = den + 2*scheme['imp_coeff'][k]*np.cos( (k+1)*kdx )
   return num / den

list_schemes = list() # List of FD schemes to plot

# ==============================================================================
#                       Explicit Standard Centred Differences
# ==============================================================================
# coefficients 2nd order explicit 3 pts
list_schemes.append( { 'exp_coeff': [2/4, -1/4]
                     , 'imp_coeff': []
                     , 'style'    : ['magenta', 'solid', '', 'STD 3pt o2']})
# coefficients 4th order explicit 5 pts
list_schemes.append( { 'exp_coeff': [6/16, -4/16, 1/16]
                     , 'imp_coeff': []
                     , 'style'    : ['green', 'solid', '', 'STD 5pt o4']})
# coefficients 6th order explicit 7 pts
list_schemes.append( { 'exp_coeff': [20/64, -15/64, 6/64, -1/64]
                     , 'imp_coeff': []
                     , 'style'    : ['orange', 'solid', '', 'STD 7pt o6']})
# coefficients 8th order explicit 9 pts
list_schemes.append( { 'exp_coeff': [70/256, -56/256, 28/256, -8/256, 1/256]
                     , 'imp_coeff': []
                     , 'style'    : ['blue', 'solid', '', 'STD 9pt o8']})
# coefficients 10th order explicit 11 pts
list_schemes.append( { 'exp_coeff': [252/1024, -210/1024, 120/1024, -45/1024, 10/1024, -1/1024]
                     , 'imp_coeff': []
                     , 'style'    : ['red', 'solid', '', 'STD 11pt o10']})
# coefficients 12th order explicit 13 pts
list_schemes.append( { 'exp_coeff': [924/4096, -792/4096, 495/4096, -220/4096, 66/4096, -12/4096, 1/4096]
                     , 'imp_coeff': []
                     , 'style'    : ['gray', 'solid', '', 'STD 13pt o12']})
# # ==============================================================================
# #                  Implicit Centred Differences (Lele JPC 1992)
# # ==============================================================================
# # 4th order family compact 3imp-5exp
# alpha=0.3; beta=0
# list_schemes.append( { 'exp_coeff': [(5+6*alpha-6*beta)/8, 0.5*(1+2*alpha+2*beta)/2, -0.5*(1-2*alpha-14*beta)/8]
#                      , 'imp_coeff': [alpha, beta]
#                      , 'style'    : ['magenta', 'dashed', '', 'CMP 3-5pt o4 a=0.3']})
# # 4th order family compact 3imp-5exp
# alpha=0.4; beta=0
# list_schemes.append( { 'exp_coeff': [(5+6*alpha-6*beta)/8, 0.5*(1+2*alpha+2*beta)/2, -0.5*(1-2*alpha-14*beta)/8]
#                      , 'imp_coeff': [alpha, beta]
#                      , 'style'    : ['green', 'dashed', '', 'CMP 3-5pt o4 a=0.4']})
# # 4th order family compact 3imp-5exp
# alpha=0.49; beta=0
# list_schemes.append( { 'exp_coeff': [(5+6*alpha-6*beta)/8, 0.5*(1+2*alpha+2*beta)/2, -0.5*(1-2*alpha-14*beta)/8]
#                      , 'imp_coeff': [alpha, beta]
#                      , 'style'    : ['orange', 'dashed', '', 'CMP 3-5pt o4 a=0.49']})
# # 6th order family compact 3imp-7exp
# alpha=0.49; beta=0
# list_schemes.append( { 'exp_coeff': [(11+10*alpha-10*beta)/16, 0.5*(15+34*alpha+30*beta)/32, 0.5*(-3+6*alpha+26*beta)/16, 0.5*(1-2*alpha+2*beta)/32]
#                      , 'imp_coeff': [alpha, beta]
#                      , 'style'    : ['blue', 'dashed', '', 'CMP 3-7pt o6 a=0.49']})
# # 8th order family compact 3imp-9exp
# alpha=0.49; beta=0
# list_schemes.append( { 'exp_coeff': [ (93+70*alpha)/128, 0.5*( 7+18*alpha)/16, 0.5*(-7+14*alpha)/32, 0.5*(1-2*alpha)/16, 0.5*(-1+2*alpha)/128]
#                      , 'imp_coeff': [alpha, beta]
#                      , 'style'    : ['gray', 'dashed', '', 'CMP 3-9pt o8 a=0.49']})
# # 10th order family compact 3imp-11exp
# alpha=0.499; beta=0
# list_schemes.append( { 'exp_coeff': [(193+126*alpha)/256,0.5*(105+302*alpha)/256,0.5*(-15+30*alpha)/64,0.5*(45-90*alpha)/512,0.5*(-5+10*alpha)/256,0.5*(1-2*alpha)/512]
#                      , 'imp_coeff': [alpha, beta]
#                      , 'style'    : ['red', 'dashed', '', 'CMP 3-11pt o10 a=0.499']})
# # ==============================================================================
# #                       Optimized Centred Differences
# # ==============================================================================
# # explicit DRP 7 pts (Bogey & Bailly 2004 JCP)
# list_schemes.append( { 'exp_coeff': [ 0.287392842460,-0.226146951809, 0.106303578770,-0.023853048191]
#                      , 'imp_coeff': []
#                      , 'style'    : ['orange', 'dotted', '', 'DRP 7pt o4']})
# # explicit4DRP 9 pts (Bogey & Bailly 2004 JCP)
# list_schemes.append( { 'exp_coeff': [0.243527493120 ,-0.204788880640 , 0.120007591680 ,-0.045211119360 , 0.008228661760]
#                      , 'imp_coeff': []
#                      , 'style'    : ['blue', 'dotted', '', 'DRP 9pt o4']})
# # explicit DRP 11 pts (Bogey & Bailly 2004 JCP)
# list_schemes.append( { 'exp_coeff': [0.234810479761700,-0.199250131285813, 0.120198310245186,-0.049303775636020, 0.012396449873964,-0.001446093078167]
#                      , 'imp_coeff': []
#                      , 'style'    : ['red', 'dotted', '', 'DRP 11pt o4']})
# # explicit DRP 13 pts (Bogey & Bailly 2004 JCP)
# list_schemes.append( { 'exp_coeff': [0.190899511506 ,-0.171503832236 , 0.123632891797 ,-0.069975429105 , 0.029662754736 ,-0.008520738659 , 0.001254597714]
#                      , 'imp_coeff': []
#                      , 'style'    : ['gray', 'dotted', '', 'DRP 13pt o4']})

#===============================================================================
#        Compute the modified wavenumber for each scheme to plot
#===============================================================================
for scheme in list_schemes:
   kstardx = transfer_function(scheme)
   if scheme['imp_coeff']:
      scheme['kstardx'] = 1-kstardx   # Transfer function
   else:
      scheme['kstardx'] = kstardx     # Transfer function

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
ax.axhline( 0, 0, np.pi, color='black')
ax.axhline( 1, 0, np.pi, color='black')

for scheme in list_schemes:
   ax.plot(kdx, scheme['kstardx'],     color=scheme['style'][0]
                                 , linestyle=scheme['style'][1]
                                 ,    marker=scheme['style'][2]
                                 ,     label=scheme['style'][3] )
# Style formatting
ax.legend(loc=2, fontsize='x-large')
ax.set_ylim([-0.1,1.1])
ax.set_xlim([0,np.pi])
ax.set_xlabel(r'$k\Delta x$', fontsize=24)
ax.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
ax.set_xticklabels(['$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$',
                     r'$\frac{3\pi}{4}$', r'$\pi$'], fontsize=24)
ax.set_ylabel(r'$D_k(k\Delta x)$', fontsize=24)
ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
ax.set_yticklabels([r'$0$', r'$0.25$', r'$0.5$', r'$0.75$', r'$1$'], fontsize=24)

plt.tight_layout()
plt.savefig('filter_transfer_function_std.pdf')

#===============================================================================
# Plot Dispersion error
# =====================
fig, ax = plt.subplots(figsize=(8,6))
ax.axhline( 1, 0, np.pi, color='black')

for scheme in list_schemes:
   ax.plot(kdx, scheme['kstardx'],     color=scheme['style'][0]
                                 , linestyle=scheme['style'][1]
                                 ,    marker=scheme['style'][2]
                                 ,     label=scheme['style'][3] )
# Style formatting
ax.legend(loc=2, fontsize='x-large')
ax.set_xlim(np.pi/16, np.pi)
ax.set_ylim(1e-5, 2)
ax.loglog()
ax.set_xlabel(r'$k\Delta x$', fontsize=24)
ax.set_xticks([np.pi/16, np.pi/8, np.pi/4, np.pi/2, np.pi])
ax.set_xticklabels([r'$\frac{\pi}{16}$', r'$\frac{\pi}{8}$', r'$\frac{\pi}{4}$',
                    r'$\frac{\pi}{2}$', r'$\pi$'], fontsize=24)
ax.set_ylabel(r'$D_k(k\Delta x)$', fontsize=24)
ax.set_yticks([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
ax.set_yticklabels([r'$10^{-6}$', r'$10^{-5}$', r'$10^{-4}$', r'$10^{-3}$',
                    r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$',], fontsize=24)
plt.tight_layout()
plt.savefig('filter_transfer_function_log_std.pdf')

# Show everything
plt.show()
