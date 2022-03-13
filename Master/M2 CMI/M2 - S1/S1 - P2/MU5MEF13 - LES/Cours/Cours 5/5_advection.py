#!/usr/bin/env python3
# finite-difference implementation for linear advection
# We are solving a_t + u a_x = 0
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import animation

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'

fig = plt.figure()
ax = plt.axes(xlim=(0, 1), ylim=(-0.2, 1.3))
sol_numer, = ax.plot([], [], lw=2, color='b')
sol_exact, = ax.plot([], [], lw=2, color='r')
time_text = ax.text(0.45, 1.175, '', fontsize=16)

plt.xlabel(r"$x$", fontsize=16)
plt.ylabel(r"$w$", fontsize=16)

def init():
   sol_numer.set_data([], [])
   sol_exact.set_data([], [])
   time_text.set_text('')
   return sol_numer, sol_exact, time_text

#===============================================================================

class FDGrid(object):

   def __init__(self, nx, ng, xmin=0.0, xmax=1.0):

      self.xmin = xmin
      self.xmax = xmax
      self.L = xmax-xmin
      self.ng = ng
      self.nx = nx

      # real data lives
      self.ilo = ng
      self.ihi = ng+nx-1

      # physical coords
      self.dx = self.L/(nx-1)
      self.x = xmin + (np.arange(nx+2*ng)-ng)*self.dx

      # storage for the solution
      self.w = np.zeros((nx+2*ng), dtype=np.float64)

      # storage for the exact solution
      self.wexact = np.zeros((nx+2*ng), dtype=np.float64)

   def scratch_array(self):
      """ return a scratch array dimensioned for our grid """
      return np.zeros((self.nx+2*self.ng), dtype=np.float64)

   def fill_BCs(self):
      """ fill a single ghost cell layer with periodic boundary conditions """
      for n in range(self.ng):
         self.w[self.ilo-(n+1)] = self.w[self.ihi-(n+1)]
         self.w[self.ihi+(n+1)] = self.w[self.ilo+(n+1)]

#===============================================================================

def solve_advection(g, a, CFL, method):
   # time info
   dt = CFL*g.dx/a
   t = 0.0

   # initialize the data
   update_exact(t)

   g.w = g.wexact.copy()

   # evolution loop
   wnew = g.scratch_array()
   if method == 'LF2': wold = g.scratch_array()
   if method == 'IM':  A = np.zeros((g.nx, g.nx), dtype=np.float64)
   # animation function.  This is called sequentially
   def animate(j):
      # fill the boundary conditions
      t=j*dt
      g.fill_BCs()

      if method == "UP": # Upwind
         for i in range(g.ilo, g.ihi+1):
            wnew[i] = g.w[i] - CFL*(g.w[i] - g.w[i-1])
      elif method == "FTCS": # Forward-Time Centred-Space
         for i in range(g.ilo, g.ihi+1):
            wnew[i] = g.w[i] - 0.5*CFL*(g.w[i+1] - g.w[i-1])
      elif method == "LW": # Lax-Wendroff
         for i in range(g.ilo, g.ihi+1):
            wnew[i] = g.w[i] - 0.5*CFL*(g.w[i+1] - g.w[i-1]) \
                             + 0.5*CFL**2*(g.w[i+1] - 2*g.w[i] + g.w[i-1])
      elif method == "BW": # Beam-Warming
         for i in range(g.ilo, g.ihi+1):
            wnew[i] = g.w[i] - 0.5*CFL*(3*g.w[i] - 4*g.w[i-1] + g.w[i-2]) \
                             + 0.5*CFL**2*(g.w[i] - 2*g.w[i-1] + g.w[i-2])
      elif method == "LF1": # Lax-Friedrichs
         for i in range(g.ilo, g.ihi+1):
            wnew[i] = 0.5*(g.w[i+1] + g.w[i-1]) - 0.5*CFL*(g.w[i+1] - g.w[i-1])
      elif method == "LF2": # Leapfrog
         if j==0:
            for i in range(g.ilo, g.ihi+1):
               wnew[i] = g.w[i] - CFL*(g.w[i] - g.w[i-1])
         else:
            for i in range(g.ilo, g.ihi+1):
               wnew[i] = wold[i] - CFL*(g.w[i+1] - g.w[i-1])
      elif method == "QUICK":
         for i in range(g.ilo, g.ihi+1):
            wnew[i] = g.w[i] - 0.125*CFL*( (3*g.w[i] + 6*g.w[i-1] - g.w[i-2]) \
                                         - (3*g.w[i+1] + 6*g.w[i] - g.w[i-1]))
      elif method == "IM":
         # create the matrix: loop over rows [ilo,ihi]
         for i in range(g.nx):
            A[i,i] = 1.0 + CFL
            A[i,i-1] = -CFL
         b = g.w[g.ilo:g.ihi+1]
         wnew[g.ilo:g.ihi+1] = np.linalg.solve(A,b)
      else:
         sys.exit("invalid method")

      # store the old solution for Leapfrog
      if method=="LF2":
         wold[:] = g.w[:]

      # store the updated solution
      g.w[:] = wnew[:]
      update_exact(t)

      sol_numer.set_data(g.x, g.w)
      sol_exact.set_data(g.x, g.wexact)

      time_text.set_text(r'$t = %s$' %(format(t,'.2f')))

      return sol_numer, sol_exact, time_text

   # call the animator.  blit=True means only re-draw the parts that have changed.
   anim = animation.FuncAnimation(fig, animate, init_func=init,
                                  frames=999999, interval=10, blit=True)
   plt.show()

#===============================================================================

case = 'tophat'   # 'tophat', 'gaussian', 'cosinehat', 'sine'
method = 'FTCS'   # 'UP', 'FTCS', 'LW', 'BW', 'LF1', 'LF2', 'QUICK', 'IM'
# create the grid
nx = 100
ng = 1
g = FDGrid(nx, ng)

# define the CFL and speed

CFL = 0.9
a = 1.0

if case == 'tophat':
   def update_exact(t):
      xx = (g.x-a*t)%g.L
      g.wexact[:] = 0
      g.wexact[np.logical_and( xx >= 1./3., xx <= 2./3.)] = 1.0
elif case == 'gaussian':
   sigma = 0.1
   def update_exact(t):
      xx = (g.x-a*t)%g.L
      g.wexact[:] = np.exp(-0.5*( (xx-0.5)/sigma)**2)
elif case == 'cosinehat':
   def update_exact(t):
      xx = (g.x-a*t)%g.L
      g.wexact[:] = [np.cos(np.pi*5/g.L*(el - g.L/10)) if el<g.L/5 else 0 for el in xx]
elif case == 'sine':
   def update_exact(t):
      xx = (g.x-a*t)%g.L
      g.wexact[:] = np.sin(2.0*np.pi*xx/g.L)

solve_advection(g, a, CFL, method)
