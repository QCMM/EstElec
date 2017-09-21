#! /usr/bin/python

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import sys

x = np.linspace(-4,4,1000)
y = np.linspace(-4,4,1000)
z = np.linspace(-4,4,1000)
X,Y = np.meshgrid(x, y)

r1 = np.linspace(0,6,1000)
r2 = np.linspace(0,6,1000)
R1,R2 = np.meshgrid(r1, r2)

zeta =  1.6875

# wave function and radial wave function
def r(x,y,z):
    return np.sqrt(x**2+y**2+z**2)

def r12(x1,y1,z1, x2, y2,z2):
    return np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

def rad_1s(r1, r2,zeta):
    return 4*(zeta**3)*np.exp(-zeta*r1-zeta*r2)

def rad_prob_1s(r1, r2, zeta):
    return rad_1s(r1,r2,zeta)**2 * r1**2*r2**2

def rad_2s(r1, r2, zeta):
    return 0.5*(zeta**3)*np.exp((-zeta*r1-zeta*r2)/2)*(1-(zeta*r1/2))*(1-(zeta*r2)/2)

def rad_prob_2s(r1, r2, zeta):
    return rad_2s(r1,r2,zeta)**2 * r1**2*r2**2

def cart_1s(x1,y1,z1, x2,y2,z2,zeta):
    return (1/np.pi)*(zeta**3)*np.exp(-zeta*(r(x1,y1,z1)+r(x2,y2,z2)))

def cart_2s(x1,y1,z1,x2,y2,z2,zeta):
    return 1/(32*np.pi)*(zeta**3)*np.exp((-zeta*r(x1,y1,z1)-zeta*r(x2,y2,z2))/2)*(2-zeta*r(x1,y1,z1))*(2-zeta*r(x2,y2,z2))

def cart_2p(x1,y1,z1,x2,y2,z2,zeta):
    return 1/(32*np.pi)*(zeta**5)*r(x1,y1,z1)*r(x2,y2,z2)*np.exp((-zeta*r(x1,y1,z1)-zeta*r(x2,y2,z2))/2)*x1/r(x1,y1,z1)*x2/r(x2,y2,z2)

def hyl1_cart(x1,y1,z1, x2,y2,z2, zeta):
    zeta = 1.81607
    zeta = zeta
    print zeta
    b = 0.364
    N = 1.295
    return (N*(1+b*r12(x1,y1,z1,x2,y2,z2))*np.exp(-zeta*(r(x1,y1,z1)+r(x2,y2,z2))))

def hyl2_cart(x1,y1,z1, x2,y2,z2, zeta):
    zetaH = 1.81607
    print zetaH
    c1 = 0.130815
    c2 = 0.29178
    N = 1.330839
    return (N*(1 + c1*(r(x1,y1,z1) - r(x2,y2,z2))**2 + c2*r12(x1,y1,z1,x2,y2,z2))*np.exp(-zetaH*(r(x1,y1,z1)+r(x2,y2,z2))))

def corr_func(x1,y1,z1, x2,y2,z2, zeta):
    return (hyl2_cart(x1,y1,z1,x2,y2,z2,zeta) - cart_1s(x1,y1,z1,x2,y2,z2,zeta))

def corr_func_1D(func,zeta):
    x1 = np.linspace(-8,8,10000)
    print x1
    y1 = 0.0
    z1 = 0.0
    x2 = 4.0
    y2 = 0.0
    z2 = 0.0
    plt.plot(x1, func(x1,y1,z1,x2,y2,z2,zeta))
    plt.xlabel("X/a_0")
    plt.ylabel("Psi_100/a_0")
    plt.show()

# Ploting Function

def plot_3D_rad(X,Y,zeta,func):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X, Y, func(X,Y,zeta), cmap=cm.coolwarm,
                              linewidth=0, antialiased=False)

    # Customize the z axis.
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

def plot_3D_cart(X,Y,z,x1,y1,z1,ze,func):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X, Y, func(X,Y,z,x1,y1,z1,ze), cmap=cm.coolwarm,
                              linewidth=0, antialiased=False)

    # Customize the z axis.
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

#plot_3D_rad(R1,R2, zeta,rad_prob_1s)
#plot_3D_rad(R1,R2,zeta,rad_prob_2s)

# Definir la grilla
''' Wavefunction of one electron in the plane when fixing the position of the second 
electron '''

x1 = np.linspace(-4,4,2000)
y1 = np.linspace(-4,4,2000)
z1 = 0.0
X1,Y1 = np.meshgrid(x1, y1)

z2 = 0.0
x2 = 0.5
y2 = 0.0

#plot_3D_cart(X1,Y1,z1,x2,y2,z2,zeta,cart_1s)
#plot_3D_cart(X1,Y1,z1,x2,y2,z2,zeta,cart_2s)
#plot_3D_cart(X1,Y1,z1,x2,y2,z2,zeta,cart_2p)

''' Hylleras wavefunction '''

#plot_3D_cart(X1,Y1,z1,x2,y2,z2,zeta,hyl1_cart)
#plot_3D_cart(X1,Y1,z1,x2,y2,z2,zeta,hyl2_cart)
#plot_3D_cart(X1,Y1,z1,x2,y2,z2,zeta,corr_func)
#plot_3D_cart(X1,Y1,z1,x2,y2,z2,zeta,hyl2_cart)
corr_func_1D(corr_func,zeta)

