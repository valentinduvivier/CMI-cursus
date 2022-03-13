#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rcParams.update({'font.size': 14})

def plotS(R, axisbox = [-4, 4, -4, 4], npts=500, scheme_name='scheme'):
   """
   Compute |R(z)| over a fine grid on the region specified by axisbox
   and do a contour plot with contourf (filled contours)
   to show the region of absolute stability.
   """
   xa, xb, ya, yb = axisbox
   x = np.linspace(xa,xb,npts)
   y = np.linspace(ya,yb,npts)
   X,Y = np.meshgrid(x,y)
   Z = X + 1j*Y
   Rval = R(Z)
   Rabs = abs(Rval)

   plt.figure()
   # plot interior, exterior, as green and white:
   levels = [-1e9,1,1e9]
   CS1 = plt.contourf(X, Y, Rabs, levels, colors = ('r', 'w'))

   # plot boundary as a black curve:
   CS2 = plt.contour(X, Y, Rabs, [1,], colors = ('k',), linewidths = (2,))

   plt.title('Absolute stability for %s' %scheme_name)
   plt.grid()
   plt.plot([xa,xb],[0,0],'k')  # x-axis
   plt.plot([0,0],[ya,yb],'k')  # y-axis
   plt.xlabel(r'$Re(\lambda \Delta t)$')
   plt.ylabel(r'$Im(\lambda \Delta t)$')
   plt.axis('scaled')  # scale x and y same so that circles are circular
   plt.axis(axisbox)   # set limits
   plt.savefig('%s.pdf'%scheme_name)
   plt.show()


# ## Theta Method
# $$w^{n+1} = w^n + k[(1-\theta)f(w^n) + \theta f(w^{n+1}]. $$
# The right hand side is a convex combination of $f(w^n)$ and $f(w^{n+1})$.  Special cases are:
#  - Forward Euler: $\theta = 0$,
#  - Trapezoid: $\theta = 1/2$,
#  - Backward Euler: $\theta = 1$.

R = lambda z: (1. + (1-theta)*z) / (1-theta*z)
theta = 0
plotS(R, scheme_name='Forward Euler')
theta = 0.5
plotS(R, scheme_name='Trapezoid')
theta = 1
plotS(R, scheme_name='Backward Euler')

R = lambda z: (1 + 5/12.*z) / (1 - 7/12.*z + 1/12.*z**2)
plotS(R, scheme_name='BDF2')
# ## Taylor series methods
# $p$th order Taylor series method gives $R(z)$ equal to the first $p+1$ terms
# of the Taylor series expansion of $e^z = \sum_{j=0}^\infty z^j/j!$
for r in range(1,10):
   def R(z):
      # return Rz = 1 + z + 0.5z^2 + ... + (1/r!) z^r
      Rz = 1.
      term = 1.
      print(r)
      for j in range(1,r+1):
         term = term * z/float(j)
         Rz = Rz + term
      return Rz
   plotS(R, npts=500,scheme_name='Taylor order r=%s'%r)
