#! /usr/bin/python

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import math

x = np.linspace(-4,4,1000)
y = np.linspace(-4,4,1000)
z = np.linspace(-4,4,1000)
X,Y = np.meshgrid(x, y)

r1 = np.linspace(0,4,1000)
r2 = np.linspace(0,4,1000)
R1,R2 = np.meshgrid(r1, r2)

Z = 2.0
zeta = Z - 5/16

# wave function and radial wave function
def phi(r1, r2):
    Z = 2.0
    zeta = Z - 5/16
    return (1/np.pi)*(zeta**3)*(np.exp(-zeta*r1)*np.exp(-zeta*r2))

def phi_rad(r1, r2):
    Z = 2.0
    zeta = Z - 5/16
    return (16/np.pi)*(zeta**3)*(np.exp(-2*zeta*r1)*np.exp(-2*zeta*r2))*r1**2*r2**2

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(R1, R2, phi_rad(R1,R2), cmap=cm.coolwarm, linewidth=0, antialiased=False)

# Customize the z axis.
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

''' Wavefunction of one electron in the plane when fixing the position of the second 
electron '''

def phi_cart(x1,y1,z1, x2,y2,z2):
    Z = 2.0
    #zeta = Z - 5/16
    zeta = 1.68
    return (1/np.pi)*(zeta**3)*(np.exp(-zeta*np.sqrt(x1**2+y1**2+z1**2))*np.exp(-zeta*np.sqrt(x2**2+y2**2+z2**2)))

z1 = 0.0
z2 = 0.0
x2 = 0.2
y2 = 0.6
x1 = np.linspace(-1,1,1000)
y1 = np.linspace(-1,1,1000)
X1,Y1 = np.meshgrid(x1, y1)

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X1, Y1, phi_cart(X1,Y1,z1,x2,y2,z2)**2, cmap=cm.coolwarm, linewidth=0, antialiased=False)

# Customize the z axis.
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

''' Hylleras wavefunction '''

def hyl_rad(r1,r2):
    ze = 1.849
    b = 0.364
    return (np.exp(-zeta*r1)*np.exp(-zeta*r2)*(1+b*(r1-r2)))


def hyl_cart(x1,y1,z1, x2,y2,z2):
    #zeta = 1.849
    zeta = 1.68
    b = 0.364
    return ((1/np.pi)*(zeta**3)*np.exp(-zeta*np.sqrt(x1**2+y1**2+z1**2))*np.exp(-zeta*np.sqrt(x2**2+y2**2+z2**2))*(1+b*np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)))

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X1, Y1, hyl_cart(X1,Y1,z1,x2,y2,z2)**2 - phi_cart(X1,Y1,z1,x2,y2,z2)**2, cmap=cm.coolwarm, linewidth=0, antialiased=False)

# Customize the z axis.
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

'''  Excited states of He '''

#1s2s

def phi_2s(x1,y1,z1,x2,y2,z2):
    z = 1.68
    return zeta**3*(2-zeta*np.sqrt(x1**2+y1**2+z1**2))*np.exp(-zeta*np.sqrt(x1**2+y1**2+z1**2))*(2-zeta*np.sqrt(x2**2+y2**2+z2**2))*np.exp(-zeta*np.sqrt(x2**2+y2**2+z2**2))

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X1, Y1, phi_2s(X1,Y1,z1,x2,y2,z2)**2, cmap=cm.coolwarm, linewidth=0, antialiased=False)

# Customize the z axis.
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()
