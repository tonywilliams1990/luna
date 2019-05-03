import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
import argparse

lines = ["-","--","-.",":"]
linecycler = cycle(lines)

data = np.loadtxt( "./DATA/Spectral_ODE_Nonlinear.dat" )
x = data[:,0]
u_exact = data[:,1]
u_spectral = data[:,2]
u_diff = data[:,3]

max_error = np.max( np.abs( u_diff ) )
print "max_error = ", max_error

plt.figure()
plt.plot( x, u_exact, 'r', label="exact" )
line = next(linecycler)
plt.plot( x, u_spectral, c='k', linestyle=line, label="spectral" )
plt.legend()

plt.figure()
plt.plot( x, u_diff, 'r--', label="difference" )
plt.legend()

plt.show()
