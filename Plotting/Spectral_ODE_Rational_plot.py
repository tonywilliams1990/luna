import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt( "./DATA/Spectral_ODE_Rational.dat" )
x = data[:,0]
u_exact = data[:,1]
u_spectral = data[:,2]
u_diff = data[:,3]

max_error = np.max( np.abs( u_diff ) )
print "max_error = ", max_error
#index = np.where( np.abs( u_diff ) == max_error )
#print "index = ", index[0][0]
#print "x location = ", x[ index[0][0] ]

plt.figure()
plt.plot( x, u_exact, 'r', label="exact" )
plt.plot( x, u_spectral, 'k', label="spectral" )
plt.legend()

axes = plt.gca()
axes.set_xlim([0,40])

#plt.figure()
#plt.plot( x, u_diff, 'r--', label="difference" )
#plt.legend()

plt.show()
