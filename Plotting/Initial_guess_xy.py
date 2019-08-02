import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt( "./DATA/Initial_guess_xy.dat" )
xy     = data[:,0]
u_g_x  = data[:,1]
u_g_y  = data[:,2]
exact  = data[:,3]


xy_max = 10
f_max = 2

plt.figure()
plt.plot( xy, u_g_x, 'k-', label="u_g_x" )
plt.plot( xy, u_g_y, 'k--', label="u_g_y" )
plt.plot( xy, exact, 'r:', label="exact" )
plt.legend()
plt.grid(True)

axes = plt.gca()
axes.set_xlim([-0.1,xy_max])
axes.set_ylim([-0.1,f_max])

plt.show()
