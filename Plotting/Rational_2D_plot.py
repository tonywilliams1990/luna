#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

N_levels = 11

data = np.loadtxt("./DATA/Rational_2D.dat")

x = data[:,0]
y = data[:,1]

u_exact = data[:,2]
u_spectral = data[:,3]
error = u_exact - u_spectral

print "max error = ", np.max( np.abs( error ) )

min_x = np.min(x)
max_x = np.max(x)
min_y = np.min(y)
max_y = np.max(y)

npts = 500

xi = np.linspace( min_x, max_x, npts )
yi = np.linspace( min_y, max_y, npts )
u_exacti = mlab.griddata( x, y, u_exact, xi, yi, interp = 'linear' )
u_spectrali = mlab.griddata( x, y, u_spectral, xi, yi, interp = 'linear' )
errori = mlab.griddata( x, y, error, xi, yi, interp = 'linear' )

origin = 'lower'
cmap = plt.cm.YlGnBu_r

levels = np.linspace( np.min( u_exact ), np.max( u_exact ), N_levels )

# Exact solution

CS = plt.contourf(xi, yi, u_exacti, levels,
                  cmap=cmap,
                  origin=origin,
                  extend='both')


CB = plt.colorbar(CS, shrink=1)
CB.set_ticks(levels)

plt.xlabel( "x" )
plt.ylabel( "y", rotation='horizontal' )
plt.title( "Exact solution" )

axes = plt.gca()
axes.set_xlim([min_x,max_x])
axes.set_ylim([min_y,max_y])

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


# Error

levels = np.linspace( np.min( error ), np.max( error ), N_levels )
plt.figure()

CS = plt.contourf(xi, yi, errori, levels,
                  cmap=cmap,
                  origin=origin,
                  extend='both')


CB = plt.colorbar(CS, shrink=1)
CB.set_ticks(levels)

plt.xlabel( "x" )
plt.ylabel( "y", rotation='horizontal' )
plt.title( "Error" )

axes = plt.gca()
axes.set_xlim([min_x,max_x])
axes.set_ylim([min_y,max_y])

plt.show()
