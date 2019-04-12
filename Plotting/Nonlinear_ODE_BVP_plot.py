import numpy as np
import matplotlib.pyplot as plt

data_exact = np.loadtxt( "./DATA/ODE_BVP_exact_solution.dat" )
data_numerical = np.loadtxt( "./DATA/ODE_BVP_numerical_solution.dat" )

x = data_exact[:,0]
f_exact = data_exact[:,1]
fdash_exact = data_exact[:,2]
f_numerical = data_numerical[:,1]
fdash_numerical = data_numerical[:,2]
f_error = np.abs( f_exact - f_numerical )
fdash_error = np.abs( fdash_exact - fdash_numerical )
max_f_error = np.max( f_error )
max_fdash_error = np.max( fdash_error )
max_error = np.max( np.concatenate( ( f_error, fdash_error ) ) )
print "max error in f = ", max_f_error
print "max error in f' = ", max_fdash_error

fig, axarr = plt.subplots(3, sharex=True, figsize=(14,8) )
fig.suptitle( "(f-1)f'' + ((2/x^2)-1)f'^2 = 1" )
axarr[0].plot( x, f_exact, 'b', label="f_exact" )
axarr[0].plot( x, f_numerical, 'b--', label="f_numerical" )
axarr[1].plot( x, fdash_exact, 'g', label="f'_exact" )
axarr[1].plot( x, fdash_numerical, 'g--', label="f'_numerical" )
axarr[2].plot( x, f_error, 'r', label="f_error" )
axarr[2].plot( x, fdash_error, 'r--', label="f'_error" )
axes = plt.gca()
axes.set_xlim( [ 0, np.sqrt(3)/2 ] )
axes.set_ylim( [ 0, max_error * 1.1 ] )
fig.legend( loc='upper left', borderaxespad=0.1 )
plt.show()
