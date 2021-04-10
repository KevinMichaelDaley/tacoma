

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import numpy as np
N=1024
def dress(yscale): 
    plt.gcf().gca().set_xticks([0,N])
    plt.gcf().gca().set_yticks([0,N])
    plt.gcf().gca().set_xticklabels(['0.4', '1.4'],fontsize=28)
    plt.gcf().gca().set_yticklabels(['0', str(yscale)],fontsize=28)
def plot_fn(X,Y,Z):
    ax=plt.gcf().gca(projection='3d')
    surf = ax.plot_surface(X, Y, Z, cmap=plt.get_cmap('magma'),
                       linewidth=1.2, antialiased=True)
    plt.gcf().gca().view_init(80,180) 

def plot_fn2(X,Y,Z):
    ax=plt.gcf().gca()
    surf = plt.imshow(Z, cmap=plt.get_cmap('magma'))
    plt.colorbar()

def plot_phase_amplitude(S, suffix, yscale,c0=0):
    amplitude=S[:,c0+2]
    phase=S[:,c0+3]
    phase=np.minimum(np.abs(phase),np.abs(1-phase))
    phase=np.abs(phase)
    AyAsBeta=amplitude.reshape([N,N])[:,:]
    PhiAsBeta=phase.reshape([N,N])[:,:]
    plt.figure()
    plot_fn2(XX,YY,Bfreq/(Ffreq/2/np.pi))
    dress(yscale)
    plt.gcf().gca().set_xlabel('$\\beta$', fontsize=55)
    plt.savefig('AyAsBeta_%s.png'%suffix)
    plt.figure()
    plot_fn2(XX,YY,PhiAsBeta/np.pi)
    dress(0)
    plt.gcf().gca().set_xlabel('$\\beta$')
    plt.gcf().gca().set_title('$\Delta\phi$', fontsize=55)
    plt.savefig("PhiAsBeta_%s.png"%suffix)


A=np.loadtxt('sweep.txt', skiprows=0)
C=np.loadtxt('eig.txt',skiprows=0)
XX,YY=np.meshgrid(np.arange(N), np.arange(N))

plot_phase_amplitude(A, "sweep", 12,c0=0)
plot_phase_amplitude(C, "sweep_linear", 12,c0=1)

#Basin=np.loadtxt('basin.txt',skiprows=0)
#phase=Basin[:,4]
#phase=np.minimum(np.abs(phase),np.abs(1-phase))
#phase=np.abs(phase)
#Z=phase<0.1
#plt.figure()
#plt.imshow(Z.reshape([50,50])[-1::-1,::1], cmap=plt.get_cmap('Greys'))
#plt.gcf().gca().set_xlabel('$u\'_{0,1}$')
#plt.gcf().gca().set_ylabel('$u\'_{0,2}$')
#plt.gcf().gca().set_xticks([0, 15])
#plt.gcf().gca().set_yticks([0, 15])
#plt.gcf().gca().set_xticklabels(['-0.5', '0.5'])
#plt.gcf().gca().set_yticklabels(['0.5', '-0.5'])
#plt.savefig('basin.pdf')

