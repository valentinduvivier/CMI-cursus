""" collection de fonctions utiles pour les notebooks
"""
#### ANIMATION
# Initialization function: plot the background of each frame
def init():
    return None
# Animation function which updates figure data.  This is called sequentially
def animate(i, line, x, y):
    line.set_data(x, y[i,:])
    return line,

def animate_double_corde(i, line, x, y):
    line.set_data(x, y[i,:])
    return line,