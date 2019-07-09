import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt( "./DATA/Spectral_Falkner.dat" )
eta  = data[:,0]
f    = data[:,1]
f_d  = data[:,2]
f_dd = data[:,3]

eta_max = 10
f_max = 2

plt.figure()
plt.plot( eta, f, 'k-', label="f" )
plt.plot( eta, f_d, 'k--', label="f'" )
plt.plot( eta, f_dd, 'k:', label="f''" )
plt.legend()

axes = plt.gca()
axes.set_xlim([0,eta_max])
axes.set_ylim([0,f_max])

plt.show()
