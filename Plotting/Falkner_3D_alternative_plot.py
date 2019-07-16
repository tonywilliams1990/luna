import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt( "./DATA/Falkner_3D_alternative.dat" )
eta  = data[:,0]
f    = data[:,1]
f_d  = data[:,2]
f_dd = data[:,3]
g    = data[:,4]
g_d  = data[:,5]
g_dd = data[:,6]

eta_max = 10
f_max = 2
g_max = 0.5


plt.figure()
plt.plot( eta, f, 'k-', label="F" )
plt.plot( eta, f_d, 'k--', label="F'" )
plt.plot( eta, f_dd, 'k:', label="F''" )
plt.legend()

axes = plt.gca()
axes.set_xlim([0,eta_max])
axes.set_ylim([0,f_max])

plt.figure()
plt.plot( eta, g, 'r-', label="G" )
plt.plot( eta, g_d, 'r--', label="G'" )
plt.plot( eta, g_dd, 'r:', label="G''" )
plt.legend()

axes = plt.gca()
axes.set_xlim([0,eta_max])
axes.set_ylim([-0.2,g_max])

plt.show()
