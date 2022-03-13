import matplotlib.pyplot as plt
import numpy.fft as fft
import numpy as np

from readvars import read_grid, read_restart, read_stats_hit

# ==================================================
#           Mise en place : ne pas toucher
# ==================================================

u0      = 10.
maille  = 5.08e-2
L_ref   = 10.8*maille/(2*np.pi)
rep = '64/'
deltat  = 3.919919198985868e-3
mu_ref  = 1.803708531754647e-5
rho_ref = 1.20720649763106
nu_ref  = mu_ref / rho_ref

# Fichier maillage
grid = rep+'grid_bl1.bin'
nx, ny, nz, x, y, z = read_grid(grid)
(xmin, xmax, ymin, ymax, zmin, zmax) = (min(x), max(x), min(y), max(y), min(z), max(z))
print("Mesh size: nx = %s, ny = %s, nz = %s"%(nx, ny, nz))

# Renormalisation après FFT
# nx*ny*nz==262144
fac = .5*L_ref / 262144**2

# Calcul k dimensionnel
kmax = int(nx/2)
k_dim = np.ndarray((kmax))
for k in range(kmax):
    k_dim[k] = (k+1) / L_ref

# Fichier x/M = 98
restart = rep+'restart_bl10056_0038.bin'
ro, rou1, rou2, rou3, roe = read_restart(restart, nx, ny, nz)

# Lecture fichier energy spectrum
es_file = rep+'es.dat'
espec = [[]]
ind1, ind2 = 1, 0
with open(es_file, "r") as f:
    for line in f.readlines():
        if ind1 > kmax:
            ind1 = 1
            ind2 += 1
            espec.append([])
        ind1 += 1
        temp = line.split()
        espec[ind2].append(float(temp[1]))

# Lecture stats.dat
file_hit = rep+'stats.dat'
dict_stats = read_stats_hit(file_hit)
for i in range(len(dict_stats['Tstar'])): dict_stats['Tstar'][i] += 42

# Comte-Bellot Corrsin : Valeurs expérimentales
# =============================================

# tU/M
tu_m_cbc = [42., 98., 171.]

# dissipation rate eps [m^2/s^3]
eps_CBC = [4740*1e-4, 633*1e-4, 174*1e-4]

# ==================================================
#               Fin de la mise en place
# ==================================================

# Question 1
L_domain  = 8.73E-2
nb_points = 64

k_max = 1 / (2 * L_domain / nb_points)
k_min = 1 / L_domain

# Question 3
# Les variables primitives pour x/M = 98 sont contenues
# dans les listes ro, rou1, rou2, rou3, roe
# Calculer à l'aide de ces listes les composantes de
# vitesse u1, u2, u3 en décommentant les lignes ci-dessous
# ==========

u1, u2, u3 = np.ndarray((nx,ny,nz)), np.ndarray((nx,ny,nz)), np.ndarray((nx,ny,nz))
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            u1[i][j][k] = rou1[i][j][k] / ro[i][j][k]
            u2[i][j][k] = rou2[i][j][k] / ro[i][j][k]
            u3[i][j][k] = rou3[i][j][k] / ro[i][j][k]

print('Finished')

# Question 4
# Calcul de la transformée de Fourier du champ de vitesse
# à l'aide de la fonction fft.fftn()
# Stocker le résultat dans u1s, u2s, u3s
# ==========

u1s = fft.fftn(u1)
u2s = fft.fftn(u2)
u3s = fft.fftn(u3)

# Question 5
# Les nombres d'ondes kx, ky, kz sont contenus dans
# les vecteurs akx, aky, akz.
#
# ak_norm est égal à la norme du vecteur k, soit
# la racine carrée de la somme des carrés de kx,
# ky et kz, que l'on peut obtenir avec akx(i)**2.
#
# n_shell est égale à l'entier le plus proche de la norme
# ak_norm.
#
# Calculer ak_norm puis n_shell dans les lignes suivantes
# Rajouter ensuite, pour k = n_shell, la somme des normes (abs())
# au carré de la vitesse dans E_spec[k].
#
# Faire ensuite la renormalisation du spectre (cf. ligne commentée)
# Décommenter les lignes suivantes
# ==========

akx, aky, akz = np.ndarray((nx)), np.ndarray((ny)), np.ndarray((nz))

for i in range(nx):
    if i<=(nx/2):
        akx[i] = i
    else:
        akx[i] = i - nx
for j in range(ny):
    if j<=(ny/2):
        aky[j] = j
    else:
        aky[j] = j - ny
for k in range(nz):
    if k<=(nz/2):
        akz[k] = k
    else:
        akz[k] = k - nz

E_spec = np.zeros((kmax))
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            ak2 = np.sqrt(akx[i]**2 + aky[i]**2 + akz[i]**2)
            n_shell = np.round(ak2)
            if n_shell>0 and n_shell<=kmax:
                n_shell = int(n_shell)
                E_spec[n_shell-1] += np.abs(u1s[i][j][k]**2 + u2s[i][j][k]**2 + u3s[i][j][k]**2)

# Renormalisation du spectre, i.e. *1/V
for k in range(kmax): E_spec[k] = E_spec[k]*fac

# E_spec is the spectrum of the kinetic energy at position x/M=98.
# We may display it :

# Question 6
# Afficher E_spec en fonction de k_dim
#
# Comparer avec les données de es.dat (contenues
# dans la liste espec[1])
#
# Décommenter le plt.show() à la fin
# Sauvegarder la figure avec plt.savefig()
# ==========

plt.figure()

plt.plot(k_dim, E_spec  , '.', color='orange', label='64 pts')
plt.plot(k_dim, espec[1], '.', color='royalblue', label='512 pts')

plt.title(r'$E_s=f(k)$ at $\frac{x}{M}=98$')
plt.xlabel(r'$k$')
plt.ylabel(r'$E_s$')

plt.xscale('log')
plt.yscale('log')

plt.grid(True, which="both")
plt.legend()

plt.tight_layout()

plt.savefig('fig 4.png')

# La comparaison des deux courbes implique de comparer 
# les énergies pour différents maillages. On remarque
# alors que pour le maillage le plus grossier les énergies
# aux faibles longeurs d'ondes ne sont pas correctement 
# retranscrites. Ainsi, la faible résolution ne permet pas
# l'étude des structures de plus petites échelles (on ne 
# les voit pas).

# De plus, même aux plus grande échelles

# Ce que je pense :
# mauvaise redistribution (ou attribution) de l'énergie :
    # étant donné le maillage grossier, les énergies normalement
    # associées aux petites échelles sont ici attribuées aux grandes échelles (k petit)
    # pour k grand, on a + d'énergie pour le maillage plus grossier, ce qui s'explique par ???

# Question 7
# Le dictionnaire dict_stats contient l'évolution de
# statistiques volumiques (<ui>,<ui^2>,...)
# dict_stats['u^2'] : variable u1^2
# dict_stats['v^2'] : variable u2^2
# dict_stats['w^2'] : variable u3^2
#
# Calculer l'évolution de l'énergie cinétique kt
# Décommenter les lignes suivantes
# ==========

Ec_t = []
for i in range(len(dict_stats['u^2'])):
    Ec_t.append(dict_stats['u^2'][i] + dict_stats['v^2'][i] + dict_stats['w^2'][i])

k_t = np.array(Ec_t) / 2

# dict_stats['Tstar'] : liste des positions tU/M = x/M
# Sauvegarder la figure avec plt.savefig()
# Décommenter les lignes suivantes
# ==========

plt.figure()

plt.plot(dict_stats['Tstar'], k_t)

plt.title('k_t as a function of adim position')
plt.xlabel(r"$\frac{tU}{M}$")
plt.ylabel(r'$<k_t>$ (m$^2$.s$^{-2}$)')

plt.grid(True)

plt.tight_layout()

plt.savefig('fig 5.png')

# It comes that energy decreases for x increasing, as seen for 512 pts case

# Question 8
# Calculer taux de dissipation epsilon
# et stocker dans epsilon_1 resultats
# deltat : pas de temps entre Ec_t[i-1] et Ec_t[i]
# Décommenter les lignes suivantes
# ==========

epsilon_1 = []
dt        = Ec_t[0] - Ec_t[1]  #L_domain / nb_points / u0
epsilon_1.append(-(k_t[1:] - k_t[:-1])/dt)

# Tracer évolution de epsilon_1 en
# fonction de dict_stats['Tstar'][1:-1]
# Décommenter les lignes suivantes
# ==========

plt.figure()

plt.plot(dict_stats['Tstar'][1:], epsilon_1[0], label='num')

plt.title('dissipation rate $\epsilon$ as a function of adim position')
plt.xlabel(r"$log(\frac{tU}{M})$")
plt.ylabel(r'$log(<\epsilon>) = - log(<\frac{d k_t}{dt})>$')

# Comparaison avec résultats expé
# Décommenter les lignes suivantes
# ==========

plt.scatter(tu_m_cbc, eps_CBC, color='r', label='exp')

plt.grid(True, which="both")
plt.legend()

plt.tight_layout()

plt.savefig('fig 6.png')

# Same evolution law, except ordonnate at the 
# origin isn't same.

# Question 9
# epsilon is dissipation rate (higher at large scales ??)
#
# it defines how much will the concentration in a place 
# dissipate across time

# Question 10
#
# D_nu(k) = 2*nu*k^2
#
# Calculer le spectre de dissipation
# à partir du spectre d'énergie E(k) et
# de la liste des k dimensionnés (k_dim)
# calculé question 6 (E_spec)
# ==========

#nu     = 1E-6
D_spec = np.zeros((kmax))
for i in range(kmax):
    D_spec[i] = 2 * nu_ref * k_dim[i]**2 * E_spec[i]

# Tracer évolution du spectre de dissipation
# Décommenter les lignes suivantes
# ==========

plt.figure()

plt.plot(k_dim, D_spec, '.', color='royalblue', label=r'$D_{\nu}(k)$')

# Ajouter spectre d'énergie dessus
# Décommenter les lignes suivantes
# ==========

plt.plot(k_dim, E_spec, '.', color='orange', label='E(k)')

plt.title(r'$D_{\nu}$ and $E$ as a function of $k$')
plt.xlabel(r"$log(k)$")
plt.ylabel(r'$log(D_{\nu})$')

plt.xscale('log')
plt.yscale('log')

plt.grid(True, which='both')
plt.legend()

plt.tight_layout()

plt.savefig('fig 7.png')

# Concordance entre dissipation de l'énergie 
# et concentration de l'énergie à différentes 
# échelles de tailles (i.e. différents k)

#PAS COMPRIS LA QUESTION + INTERPRETATION RESTE A FAIRE

# Question 11
# Calcul de l'enstrophie en fonction du temps et
# stockage dans enstrophie.
#
# Les 3 composantes de la vorticité au carré moyenné
# sont contenues dans la liste dict_stats['b2^1'],
# dict_stats['b3^2'] et dict_stats['b3^2'][i]
# Décommenter les lignes suivantes
# ==========

enstrophie = []
for i in range(len(dict_stats['b1^2'])):
    enstrophie.append((dict_stats['b1^2'][i] + dict_stats['b2^2'][i] + dict_stats['b3^2'][i]) / 2)

# Calcul de epsilon et stockage dans
# epsilon_2. La viscosité cinématique
# nu est contenue dans la variable nu_ref
# Décommenter les lignes suivantes
# ==========

epsilon_2 = []
for i in range(len(enstrophie)):
    epsilon_2.append(2 * nu_ref * enstrophie[i])

# Question 12
# Tracer l'évolution de epsilon_2 en fonction
# de dict_stats['Tstar']
#
# Comparer avec epsilon_1 et les résultats
# expérimentaux (cf. question 8)
# Décommenter les lignes suivantes
# ==========

plt.figure()

plt.plot(dict_stats['Tstar']    , epsilon_2   , '.', color='royalblue', label='epsilon_2')
plt.plot(dict_stats['Tstar'][1:], epsilon_1[0], '.', color='orange', label='epsilon_1')

plt.scatter(tu_m_cbc, eps_CBC, color='r', label='CBC')

plt.title(r'Dissipation rate for two approaches along flow')
plt.xlabel(r"$log(\frac{tU}{M})$")
plt.ylabel(r'$log(<\epsilon>) = log(\nu <\omega_i^2>)$')

plt.xscale('log')
plt.yscale('log')

plt.grid(True, which='both')
plt.legend()

plt.tight_layout()

plt.savefig('fig 8.png')

# As mentioned, dissipation not well reproduced. Then, for 64 points,
# we lack precision or date to fit theory well (num here).

# Again, same slope but not same ordonate at the origin.
# It may be we have neglected some terms (e.g. convective)
# or that we in fact are influenced by laminar flow or else
# !!!!!