import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt( "./DATA/ODE_BVP_System.dat" )
x       = data[:,0]
u       = data[:,1]
v       = data[:,2]
u_exact = data[:,3]
v_exact = data[:,4]

x_max = 10
x_max_error = 100
u_max = 2

plt.figure()
plt.plot( x, u, 'k-', label="u" )
plt.plot( x, v, 'k--', label="v" )
plt.plot( x, u_exact, 'b-', label="u_exact" )
plt.plot( x, v_exact, 'b--', label="v_exact" )
plt.legend()

axes = plt.gca()
axes.set_xlim([0,x_max])
axes.set_ylim([0,u_max])

plt.figure()
plt.plot( x, u - u_exact, 'r-', label="u_error" )
plt.plot( x, v - v_exact, 'r--', label="v_error" )
plt.legend()

axes = plt.gca()
axes.set_xlim([0,x_max_error])

plt.show()
