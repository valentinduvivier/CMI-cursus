import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

def lorenz_deriv( state, t0, sigma=10., beta=8./3, rho=35.0):
   """Compute the time-derivative of a Lorenz system."""
   x, y, z = state
   return [sigma * (y - x), x * (rho - z) - y, x * y - beta * z]

# Three trajectories
N_trajectories = 4
eps = 1e-5
x0 = [ [0.5    , 0.1    , 0.3    ],
       [0.5+eps, 0.1    , 0.3    ],
       [0.5    , 0.1+eps, 0.3    ],
       [0.5    , 0.1    , 0.3+eps]]

colors=['forestgreen','b','yellow','r']

# Solve for the trajectories
T = 100
N = 10000
t = np.linspace(0, T, N)
x_t = np.asarray([integrate.odeint(lorenz_deriv, x0i, t) for x0i in x0])

# Set up figure & 3D axis for animation
fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection='3d')

# set up lines and points
lines = sum([ax.plot([], [], [], '-', c=c) for c in colors], [])
pts   = sum([ax.plot([], [], [], 'o', c=c,markeredgecolor="k") for c in colors], [])

# prepare the axes limits
ax.set_xlim((-25, 25))
ax.set_ylim((-35, 35))
ax.set_zlim((  5, 55))

# set point-of-view: specified by (altitude degrees, azimuth degrees)
ax.view_init(30, 135)

# initialization function: plot the background of each frame
def init():
   for line, pt in zip(lines, pts):
      line.set_data([], [])
      # line.set_3d_properties([])

      pt.set_data([], [])
      pt.set_3d_properties([])
   return lines + pts

# animation function. This will be called sequentially with the frame number
def animate(i):
    for line, pt, xi in zip(lines, pts, x_t):
        x, y, z = xi[:i].T
        line.set_data(x, y)
        line.set_3d_properties(z)

        pt.set_data(x[-1:], y[-1:])
        pt.set_3d_properties(z[-1:])

    # ax.view_init(30, 0.3 * i)
    # fig.canvas.draw()
    return lines + pts

# instantiate the animator.
anim = animation.FuncAnimation(fig, animate, frames=10000,
                               interval=2, blit=True)

# anim.save('lorentz_attractor.mp4', fps=15, extra_args=['-vcodec', 'libx264'])
plt.show()
