import numpy as np
import pylab as plt

from ipywidgets import interact



# Frequence d'echantillonnage
Fe = 500
# vecteur frequence
f = np.linspace(0,Fe,1000)
# Construction du cercle unite
z = np.exp(-2*1j*np.pi*f/Fe)


def RIF_1(z1=-1):
    global z
    H = (z-z1)

    # Representation graphique
    plt.figure(figsize=(7,4))
    plt.subplot(121)
    plt.plot(np.real(z),np.imag(z),linewidth=2)
    plt.plot(np.real(z1),np.imag(z1),'or')
    plt.title(u"Position du zéro dans le plan complexe")
    plt.grid()
    plt.axis('equal')
    plt.axis([-1.5,1.5,-1.1,1.1])

    #plt.figure(figsize=(6,6))
    plt.subplot(122)
    plt.plot(f/Fe,abs(H),lw=3)
    plt.xlabel(r"$f/F_E$")
    plt.ylabel(r"$|H(e^{2j\pi f T_e})|$")
    plt.title(u"Gain en fréquences")
    plt.axis([0,1,0,3])
    plt.grid()

def RIF2(z1m, theta1):
    z1 = z1m * np.exp(1j*theta1)
    z2 = z1m * np.exp(-1j*theta1)
    
    H=(z-z1)*(z-z2)

    plt.figure(figsize=(7,4))
    plt.subplot(121)
    plt.plot(np.real(z),np.imag(z),linewidth=2)
    plt.plot(np.real(z1),np.imag(z1),'or')
    plt.plot([0, np.real(z1)],[0, np.imag(z1)], lw=1, color='k')
    plt.plot(np.real(z2),np.imag(z2),'ob')
    plt.plot([0, np.real(z2)],[0, np.imag(z2)], lw=1, color='k')
    plt.title(u"Position des zéros dans le plan complexe")
    plt.grid()
    plt.axis('equal')
    plt.axis([-1.1,1.1,-1.1,1.1])

    # Representation graphique
    plt.subplot(122)
    plt.plot(f/Fe,abs(H),lw=3)
    plt.xlabel(r"$f/F_E$")
    plt.ylabel(r"$|H(e^{2j\pi f T_e})|$")
    plt.title(u"Gain en fréquences")
    plt.axis([0,1,0,5])
    plt.grid()


def RII_1(p1=-1,zz_arg=0):
    
    # Fonction de transfert en z
    H = 1/(z-p1)
    zz = np.exp(1j*zz_arg)

    # Representation graphique
    plt.figure(figsize=(7,4))
    plt.subplot(121)
    plt.plot(np.real(z),np.imag(z),linewidth=2)
    plt.plot(np.real(p1),np.imag(p1),'or')
    #plt.plot([0, np.real(p1)],[0, np.imag(p1)], lw=1, color='k')
    plt.plot([np.real(zz), np.real(p1)],[np.imag(zz), np.imag(p1)], lw=1, color='k')
    plt.plot(np.real(zz),np.imag(zz),'ok')
    
    plt.title(u"Position du zéro dans le plan complexe")
    plt.grid()
    plt.axis('equal')
    plt.axis([-1.1,1.1,-1.1,1.1])

    #plt.figure(figsize=(6,6))
    plt.subplot(122)
    plt.plot(f/Fe,abs(H),lw=3)
    plt.plot(zz_arg/(2*np.pi)*np.ones_like(f), np.linspace(-np.max(np.abs(H)),np.max(np.abs(H)), len(f)))
    plt.xlabel(r"$f/F_E$")
    plt.ylabel(r"$|H(e^{2j\pi f T_e})|$")
    plt.title(u"Gain en fréquences")
    plt.axis([0,1.0,0,np.max(np.abs(H))])
    plt.grid()

def RII_2_rec(p1_mod=-1,p1_arg=-1,zz_arg=0):
    
    p1 = p1_mod * np.exp(1j * p1_arg)
    p2 = p1_mod * np.exp(-1j * p1_arg)
    # Fonction de transfert en z
    H = 1/((z-p1)*(z-p2))
    zz = np.exp(1j*zz_arg)

    # Representation graphique
    plt.figure(figsize=(7,4))
    plt.subplot(121)
    plt.plot(np.real(z),np.imag(z),linewidth=2)
    plt.plot(np.real(p1),np.imag(p1),'or')
    plt.plot(np.real(p2),np.imag(p2),'or')
    #plt.plot([0, np.real(p1)],[0, np.imag(p1)], lw=1, color='k')
    plt.plot([np.real(zz), np.real(p1)],[np.imag(zz), np.imag(p1)], lw=1, color='k')
    plt.plot([np.real(zz), np.real(p2)],[np.imag(zz), np.imag(p2)], lw=1, color='k')
    plt.plot(np.real(zz),np.imag(zz),'ok')
    
    plt.title(u"Position du zéro dans le plan complexe")
    plt.grid()
    plt.axis('equal')
    plt.axis([-1.1,1.1,-1.1,1.1])

    #plt.figure(figsize=(6,6))
    plt.subplot(122)
    plt.plot(f/Fe,abs(H),lw=3)
    plt.plot(zz_arg/(2*np.pi)*np.ones_like(f), np.linspace(-np.max(np.abs(H)),np.max(np.abs(H)), len(f)))
    plt.xlabel(r"$f/F_E$")
    plt.ylabel(r"$|H(e^{2j\pi f T_e})|$")
    plt.title(u"Gain en fréquences")
    plt.axis([0,1.0,0,np.max(np.abs(H))])
    plt.grid()

def RII_2(p1_mod=-1,p1_arg=-1, z1_mod=-1,z1_arg=-1,zz_arg=0):
    
    p1 = p1_mod * np.exp(1j * p1_arg)
    p2 = p1_mod * np.exp(-1j * p1_arg)
    z1 = z1_mod * np.exp(1j * z1_arg)
    z2 = z1_mod * np.exp(-1j * z1_arg)
    
    # Fonction de transfert en z
    H = ((z-z1)*(z-z2))/((z-p1)*(z-p2))
    zz = np.exp(1j*zz_arg)

    # Representation graphique
    plt.figure(figsize=(7,4))
    plt.subplot(121)
    plt.plot(np.real(z),np.imag(z),linewidth=2)
    plt.plot(np.real(p1),np.imag(p1),'or')
    plt.plot(np.real(p2),np.imag(p2),'or')
    plt.plot(np.real(z1),np.imag(z1),'ob')
    plt.plot(np.real(z2),np.imag(z2),'ob')
    #plt.plot([0, np.real(p1)],[0, np.imag(p1)], lw=1, color='k')
    plt.plot([np.real(zz), np.real(p1)],[np.imag(zz), np.imag(p1)], lw=1, color='k')
    plt.plot([np.real(zz), np.real(p2)],[np.imag(zz), np.imag(p2)], lw=1, color='k')
    plt.plot([np.real(zz), np.real(z1)],[np.imag(zz), np.imag(z1)], lw=1, color='r')
    plt.plot([np.real(zz), np.real(z2)],[np.imag(zz), np.imag(z2)], lw=1, color='r')
    plt.plot(np.real(zz),np.imag(zz),'ok')
    plt.axis('equal')
    plt.title(u"Position du zéro dans le plan complexe")
    plt.grid()
    plt.axis([-1.1,1.1,-1.1,1.1])
    #plt.savefig("Positions.pdf")
    #plt.figure(figsize=(6,6))
    plt.subplot(122)
    plt.plot(f/Fe,abs(H),lw=3)
    plt.plot(zz_arg/(2*np.pi)*np.ones_like(f), np.linspace(-np.max(np.abs(H)),np.max(np.abs(H)), len(f)))
    plt.xlabel(r"$f/F_E$")
    plt.ylabel(r"$|H(e^{2j\pi f T_e})|$")
    plt.title(u"Gain en fréquences")
    plt.axis([0,1.0,0,np.max(np.abs(H))])
    plt.grid()
    #plt.savefig("Gain.pdf")

