#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# plt.rc('xtick', labelsize=14)
# plt.rc('ytick', labelsize=14)
# plt.rc('font', size=22)
plt.rc('text', usetex=True)
plt.rcParams.update({'font.size': 14})
def exact(N,R):
   return [(np.exp(j*R)-1)/(np.exp(N*R)-1) for j in range(N+1)]

def central(N,R):
   term = (2+R)/(2-R)
   return [(term**j-1)/(term**N-1) for j in range(N+1)]

def upwind(N,R):
   term = 1+R
   return [(1 - term**j)/(1 - term**N) for j in range(N+1)]

N=100
j = np.linspace(0,1,N+1)
plt.figure(figsize=(6,3))
# R = R_L/N
plt.plot(j, exact(N,  1/N), color='k'      , ls='-', label=r'$R_m=0.01, R_L=1$'  )
plt.plot(j, exact(N,  5/N), color='orange' , ls='-', label=r'$R_m=0.05, R_L=5$'  )
plt.plot(j, exact(N, 10/N), color='b'      , ls='-', label=r'$R_m=0.10, R_L=10$' )
plt.plot(j, exact(N, 50/N), color='magenta', ls='-', label=r'$R_m=0.50, R_L=50$' )
plt.plot(j, exact(N,100/N), color='green'  , ls='-', label=r'$R_m=1.00, R_L=100$')
plt.plot(j, exact(N,200/N), color='r'      , ls='-', label=r'$R_m=2.00, R_L=200$')

plt.plot(j[-2], exact(N,  1/N)[-2], color='k'      , marker='o')
plt.plot(j[-2], exact(N,  5/N)[-2], color='orange' , marker='o')
plt.plot(j[-2], exact(N, 10/N)[-2], color='b'      , marker='o')
plt.plot(j[-2], exact(N, 50/N)[-2], color='magenta', marker='o')
plt.plot(j[-2], exact(N,100/N)[-2], color='green'  , marker='o')
plt.plot(j[-2], exact(N,200/N)[-2], color='r'      , marker='o')

plt.title(r'Exact solution, $N=%s$'%N)
plt.xlim([0,1.1])
plt.ylim([-0.05,1])
plt.grid()
plt.legend()
plt.show()
#==================================================================
# Upwind, fine mesh vs coarse mesh
R = 0.25
N = 60
j = np.linspace(0,1,N+1)
plt.figure(figsize=(4,4))
plt.plot(j,  exact(N,R), color='k', ls='-', label='exact')
plt.plot(j, upwind(N,R), color='b', ls='-', label='upwind (fine mesh)')
N = int(N/4)
j = np.linspace(0,1,N+1)
plt.plot(j, upwind(N,4*R), color='b', ls='--', label='upwind (coarse mesh)')
plt.xlim([    0,1])
plt.ylim([-0.15,1])
plt.grid()
plt.ylabel(r'$y/\delta$')
plt.xlabel(r'$U/U_e$')
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig('bl_mesh_upwind.pdf')

#==================================================================
# Central, fine mesh vs coarse mesh
R = 0.25
N = 60
j = np.linspace(0,1,N+1)
plt.figure(figsize=(4,4))
plt.plot(j,  exact(N,R), color='k', ls='-', label='exact')
plt.plot(j, central(N,R), color='r', ls='-', label='central (fine mesh)')
N = int(N/4)
j = np.linspace(0,1,N+1)
plt.plot(j, central(N,4*R), color='r', ls='--', label='central (coarse mesh)')
plt.xlim([    0,1])
plt.ylim([-0.15,1])
plt.grid()
plt.ylabel(r'$y/\delta$')
plt.xlabel(r'$U/U_e$')
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig('bl_mesh_central.pdf')

#==================================================================
# Comparison Central vs Upwind coarse mesh
R = 0.25
N = 60
j = np.linspace(0,1,N+1)
plt.figure(figsize=(4,4))
plt.plot(j,  exact(N,R), color='k', ls='-', label='exact')
N = int(N/4)
j = np.linspace(0,1,N+1)
plt.plot(j,  upwind(N,4*R), color='b', ls='--', label='upwind (coarse mesh)')
plt.plot(j, central(N,4*R), color='r', ls='--', label='central (coarse mesh)')
plt.xlim([    0,1])
plt.ylim([-0.15,1])
plt.grid()
plt.ylabel(r'$y/\delta$')
plt.xlabel(r'$U/U_e$')
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig('bl_mesh_comparison.pdf')

#==================================================================
# Upwind vs central
N = 10
j = np.linspace(0,1,N+1)
R = 2.5
plt.figure(figsize=(4,3))
plt.plot(j,   exact(N,R), color='k', ls='-', label='exact')
plt.plot(j,  upwind(N,R), color='r', ls='-', label='upwind')
plt.plot(j, central(N,R), color='b', ls='-', label='central')
plt.title(r'$Re_m = %s$' %R)
plt.xlim([    0,1])
plt.ylim([-0.15,1])
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('Re_cell%s.pdf'%R)

plt.show()
