import matplotlib.pyplot as plt
import numpy as np

# ==================================================
#           Mise en place : ne pas toucher
# ==================================================

# Comte-Bellot Corrsin : Valeurs expérimentales
# =============================================

# Vitesse moyenne [m/s]
u0 = 10.0

# tU/M = x/M
tu_m_cbc = [42., 98., 171.]

# urms = sqrt(<u1^2>) [cm/s]
urms_CBC = [22.2, 12.8, 8.95]

# U^2 / u'^2
U2_u2_cbc = []
for elmt in urms_CBC: U2_u2_cbc.append(u0**2/elmt**2*1e4)

# Tu (%)
Tu_cbc = []
for elmt in U2_u2_cbc: Tu_cbc.append((elmt)**(-0.5)*100)

# Longueur intégrale transversale [m]
Lf_CBC = [2.4, 3.45, 4.9]

# Spectres d'énergie 3D :
# nombre d'onde k [cm^-1]
kk1  = [0.20,0.25,0.30,0.40,0.50,0.70,1.00,1.50,2.00,2.50,3.00,4.00,6.00,8.00,10.00,12.50,15.00,17.50,20.00]
kk2  = [0.20,0.25,0.30,0.40,0.50,0.70,1.00,1.50,2.00,2.50,3.00,4.00,6.00,8.00,10.00,12.50,15.00,17.50,20.00]
kk3  = [0.15,0.20,0.25,0.30,0.40,0.50,0.70,1.00,1.50,2.00,2.50,3.00,4.00,6.00,8.00,10.00,12.50,15.00]

# E(k,t) [cm^3/sec^2] à tUo/M=42
Ek1 = [129.,230.,322.,435.,457.,380.,270.,168.,120,89.0,70.3,47.0,24.7,12.6,7.42,3.96,2.33,1.34,0.80]

# E(k,t) [cm^3/sec^2] à tUo/M=98
Ek2 = [106.,196.,195.,202.,168.,127.,79.2,47.8,34.6,28.6,23.1,14.3,5.95,2.23,0.900,0.363,0.162,0.0660,0.0330]

# E(k,t) [cm^3/sec^2] à tUo/M=171
Ek3 = [49.7,92.0,120,125,98.0,81.5,60.2,39.4,24.1,16.5,12.5,9.12,5.62,1.69,0.520,0.161,0.0520,0.0141]
for i in range(len(kk1)):
    kk1[i] = kk1[i]*1e2
    kk2[i] = kk2[i]*1e2
    Ek1[i] = Ek1[i]*1e-6
    Ek2[i] = Ek2[i]*1e-6
for i in range(len(kk3)):
    kk3[i] = kk3[i]*1e2
    Ek3[i] = Ek3[i]*1e-6

# Simulation numérique 512^3
# ==========================
rep = '512/'

# Statistiques
res = {'Tstar': [], 'urms': []}
L_f = []
file = rep + "stat1.dat"
with open(file, "r") as f:
    for line in f.readlines():
        save = line.split()
        res['Tstar'].append(float(save[0]))
        res['urms'].append(float(save[4]))
        L_f.append(float(save[5])*1e2)

# Lecture fichier energy spectrum
es_file = rep+'es.dat'
kmax = 512//2
espec,k = [[]], [[]]
ind1, ind2 = 1, 0

with open(es_file, "r") as f:
    for line in f.readlines():
        if ind1 > kmax:
            ind1 = 1
            ind2 += 1
            k.append([])
            espec.append([])
        ind1 += 1
        temp = line.split()
        k[ind2].append(float(temp[0]))
        espec[ind2].append(float(temp[1]))

Tu, tU_m = [], []
for i in range(len(res['Tstar'])):
    Tu.append((res['urms'][i])**0.5/u0*100)
    tU_m.append(res['Tstar'][i]+42)

# ==================================================
#               Fin de la mise en place
# ==================================================


# ==================================================
#            Partie réponse aux questions
# ==================================================


# Question 1
# tU_m : liste avec positions x/M = t*u0/M
# Tu : liste avec l'intensité turbulente aux positions
# Compléter le plt.plot() en ajoutant les données à tracer :
#     tU_m pour x et Tu pour y
# Décommenter le code dessous pour tracer évolution Tu
# Décommenter le plt.show() à la fin
# ==========

def func(x, a, b, c):
    
    """ Function returning num. fit for y=a*x curve """
    
    return b * (x-c)**a

# --------------------------
# Deduce resistance value
# --------------------------

import scipy.optimize as optimization
a, b, c = optimization.curve_fit(func, np.array(tU_m), Tu)[0]
#a = -nd/2; #b = 1/sqrt(A); #c=x0/M

x  = np.log(np.linspace(tU_m[0], tU_m[-1], 100))

# ---
# coeff analytical law

n_d = -2*a
A   = 1/b**2
x0  = c      #x0==x0/M

# ------------------

plt.figure()

plt.plot(tU_m, Tu, '.', color='orange', label='exp values')
plt.plot(np.array(tU_m), b*(np.array(tU_m)-c)**a, '-', color='royalblue', label='num interp')

plt.xlabel(r"$log(\frac{tU}{M})$")
plt.ylabel(r"$log(T_u)$")
plt.title(r'Evolution de $T_u$ en fonction de  $\frac{tU}{M}$')

plt.xscale('log')
plt.yscale('log')

plt.grid(True, which="both")
plt.legend()

# # Sauvegarder figure pour rapport
# # ==========

plt.tight_layout()

# plt.savefig('fig 1b.png')

## ==========

# Question 5
# Décommenter les lignes suivantes pour tracer les spectres
# ==========

plt.figure()
for i in range(0, ind2+1):
    if np.round(res['Tstar'][i]+42,0)==42 or np.round(res['Tstar'][i]+42,0)==98 or np.round(res['Tstar'][i]+42,0)==171:
        plt.plot(k[i][:], espec[i][:], '-', label=r'$\frac{x}{M}=$%i'%(np.round(res['Tstar'][i]+42,0)))

plt.xscale('log')
plt.yscale('log')

plt.title('Spectres num vs ana')
plt.ylabel(r"$log(E_{spec})$")
plt.xlabel(r"$log(k)$)")

# ---------

L_domain  = 8.73E-2
nb_points = 512

k_max = 1 / ((2 * L_domain) / nb_points)
k_min = 1 / L_domain

# Question 7
# Les spectres expérimentaux sont contenus dans les listes Ek1, Ek2 et Ek3,
# avec kk1, kk2, kk3 associé pour les nombres d'ondes
# Utiliser plt.scatter(kk1,Ek1) pour afficher juste les points  sur la
# figure précédente (fonctionne comme plt.plot())
# ==========

plt.scatter(kk1,Ek1)
plt.scatter(kk2,Ek2)
plt.scatter(kk3,Ek3)

plt.grid(True, which="both")
plt.legend()

plt.tight_layout()

# plt.savefig('fig 2.png')

# Question 8
# tU_m : liste avec positions x/M = t*u0/M
# L_f  : liste avec l'échelle intégrale pour chaque position
#
# Tracer l'évolution de l'échelle intégrale en fonction de
# la position en aval de la grille, namely tU_m
#
# Les listes tu_m_cbc et Lf_CBC contiennent les données expé.
# A ajouter avec plt.scatter()
# Décommenter les lignes suivantes
# ==========

plt.figure()

plt.plot(tU_m, L_f, '.-', label='num')
plt.plot(tu_m_cbc, Lf_CBC, '.', label='ana')

plt.title('Echelle integrale en fonction de la position adim')
plt.xlabel(r"$\frac{tU}{M}$")
plt.ylabel(r"$L_f$")

plt.grid(True)
plt.legend(loc=2)

plt.tight_layout()

# plt.savefig('fig 3.png')