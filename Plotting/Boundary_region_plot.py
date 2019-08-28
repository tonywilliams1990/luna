#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

#TODO input beta as a parameter and make nice titles etc

N_levels = 11

data = np.loadtxt("./DATA/Boundary_region.dat")

x = data[:,0]
y = data[:,1]

u_spectral = data[:,2]

min_x = np.min(x)
max_x = np.max(x)
min_y = np.min(y)
max_y = np.max(y)

npts = 500

xi = np.linspace( min_x, max_x, npts )
yi = np.linspace( min_y, max_y, npts )
u_spectrali = mlab.griddata( x, y, u_spectral, xi, yi, interp = 'linear' )

origin = 'lower'
cmap = plt.cm.YlGnBu_r

# Spectral solution

levels = np.linspace( np.min( u_spectral ), np.max( u_spectral ), N_levels )
plt.figure()

CS = plt.contourf(xi, yi, u_spectrali, levels,
                  cmap=cmap,
                  origin=origin,
                  extend='both')


CB = plt.colorbar(CS, shrink=1)
CB.set_ticks(levels)

plt.xlabel( "x" )
plt.ylabel( "y", rotation='horizontal' )
plt.title( "Spectral solution" )

axes = plt.gca()
axes.set_xlim([min_x,max_x])
axes.set_ylim([min_y,max_y])

# Base solution (Falkner-Skan)
'''
data    = np.loadtxt( "./DATA/Boundary_region_base.dat" )
eta     = data[:,0]
UB      = data[:,1]
UBd     = data[:,2]
PhiB    = data[:,3]
ThetaB  = data[:,4]
ThetaBd = data[:,5]
PsiB    = data[:,6]

eta_max = 10
f_max = 2

plt.figure()
plt.plot( eta, UB, 'k-', label="UB" )
plt.plot( eta, UBd, 'k--', label="UB'" )
plt.plot( eta, PhiB, 'k:', label="PhiB" )
plt.plot( eta, ThetaB, 'r-', label="ThetaB" )
plt.plot( eta, ThetaBd, 'r--', label="ThetaBd" )
plt.plot( eta, PsiB, 'b-', label="PsiB" )
plt.legend()

axes = plt.gca()
axes.set_xlim([0,eta_max])
axes.set_ylim([0,f_max])
'''

plt.show()
