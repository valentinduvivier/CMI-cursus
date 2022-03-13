# Script to plot second derivative modified wavenumber
import numpy as np
import matplotlib.pyplot as plt

# main program
npt = 200                           # number of points
kdx = np.linspace( 0., np.pi, npt )

list_schemes = list()

list_schemes.append( { '(kstardx)^2': -2*(np.cos(kdx)-1)
                     , 'style'    : ['red', 'solid', '', r'$w_{j+1}-2w_j+w_{j-1}$']})
list_schemes.append( { '(kstardx)^2': np.sin(kdx)**2
                     , 'style'    : ['blue', 'solid', '', r'$(w_{j+2}-2w_j+w_{j-2})/4$']})
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
ax.plot( kdx, (kdx)**2, color='black', ls='solid', label='Exact' )  # Exact solution
for scheme in list_schemes:
   ax.plot(kdx, scheme['(kstardx)^2'],     color=scheme['style'][0]
                                     , linestyle=scheme['style'][1]
                                     ,    marker=scheme['style'][2]
                                     ,     label=scheme['style'][3] )

# Style formatting
ax.legend(loc=2, fontsize='x-large')
ax.set_ylim([0,np.pi**2])
ax.set_xlim([0,np.pi])
ax.set_xlabel(r'$k\Delta x$', fontsize=24)
ax.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
ax.set_xticklabels(['$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$',
                     r'$\frac{3\pi}{4}$', r'$\pi$'], fontsize=24)
ax.set_ylabel(r'$(k^*\Delta x)^2$', fontsize=24)
ax.set_yticks([0, np.pi**2/2, np.pi**2])
ax.set_yticklabels([r'$0$', r'$\frac{\pi^2}{2}$', r'$\pi^2$'], fontsize=24)

plt.tight_layout()
plt.savefig('second_der_modified_wavenumber.pdf')

#===============================================================================
# Plot phase error
# ================
fig, ax = plt.subplots(figsize=(8,6))
ax.plot( kdx, (kdx)**2/(kdx)**2, color='black', label='Exact' )  # Exact solution
for scheme in list_schemes:
   ax.plot(kdx, scheme['(kstardx)^2']/(kdx)**2,     color=scheme['style'][0]
                                              , linestyle=scheme['style'][1]
                                              ,    marker=scheme['style'][2]
                                              ,     label=scheme['style'][3] )

ax.legend(loc=3, fontsize='x-large')
ax.set_ylim([0,1.05])
ax.set_xlim([0,np.pi])
ax.set_xlabel(r'$k\Delta x$', fontsize=24)
ax.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
ax.set_xticklabels(['$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$',
                     r'$\frac{3\pi}{4}$', r'$\pi$'], fontsize=24)
ax.set_ylabel(r'$(k^*\Delta x)^2/ (k\Delta x)^2$', fontsize=24)

plt.tight_layout()
plt.savefig('second_der_amp_error.pdf')

# Show everything
plt.show()