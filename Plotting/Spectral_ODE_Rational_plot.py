import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt( "./DATA/Spectral_ODE_Rational.dat" )
y = data[:,0]
u_exact = data[:,1]
u_spectral = data[:,2]
u_diff = data[:,3]

#y_max = np.max( y )
y_max = 20

max_error = np.max( np.abs( u_diff ) )
print "max_error = ", max_error
index = np.where( np.abs( u_diff ) == max_error )
print "y location of max_error = ", y[ index[0][0] ]

plt.figure()
plt.plot( y, u_exact, 'r', label="exact" )
plt.plot( y, u_spectral, 'k', label="spectral" )
plt.legend()

axes = plt.gca()
axes.set_xlim([0,y_max])

plt.figure()
plt.plot( y, u_diff, 'r--', label="difference" )
plt.legend()

axes = plt.gca()
axes.set_xlim([0,y_max])

plt.show()
