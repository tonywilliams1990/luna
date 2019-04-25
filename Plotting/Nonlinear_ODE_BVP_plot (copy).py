import numpy as np
import matplotlib.pyplot as plt

#n = 4

ns = [4, 6, 8, 10]

plt.figure()

for n in ns:
    data = np.loadtxt( "./DATA/Spectral_ODE_BVP_n_" + str( n ) + ".dat" )

    x = data[:,0]
    u_exact = data[:,1]
    u_spectral = data[:,2]

    plt.plot( x, u_exact, 'r' )
    plt.plot( x, u_spectral, 'b*')

plt.show()
