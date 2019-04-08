import numpy as np
import matplotlib.pyplot as plt

data_exact = np.loadtxt( "./DATA/ODE_BVP_exact_solution.dat" )
data_numerical = np.loadtxt( "./DATA/ODE_BVP_numerical_solution.dat" )

x = data_exact[:,0]
y_exact = data_exact[:,1]
ydash_exact = data_exact[:,2]
y_numerical = data_numerical[:,1]
ydash_numerical = data_numerical[:,2]
y_error = np.abs( y_exact - y_numerical )
ydash_error = np.abs( ydash_exact - ydash_numerical )
max_y_error = np.max( y_error )
max_ydash_error = np.max( ydash_error )
max_error = np.max( np.concatenate( ( y_error, ydash_error ) ) )
print "max_error = ", max_error

f, axarr = plt.subplots(3, sharex=True)
f.suptitle( "y'' + 4y = 0" )
axarr[0].plot( x, y_exact, 'b' )
axarr[0].plot( x, ydash_exact, 'b--' )
axarr[1].plot( x, y_numerical, 'b' )
axarr[1].plot( x, ydash_numerical, 'b--' )
axarr[2].plot( x, y_error, 'r' )
axarr[2].plot( x, ydash_error, 'r--' )
axes = plt.gca()
axes.set_xlim( [ 0, np.pi / 4 ] )
axes.set_ylim( [ 0, max_error * 1.1 ] )
plt.show()
