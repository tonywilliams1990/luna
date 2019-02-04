#!/usr/bin/env python
import os
import subprocess
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import cycle

slowdown = 5
J = 50
N = 500
data = np.loadtxt( "./Heat_equation_J_"+ str(J) + "_N_" + str(N) + ".dat" )

x = data[:,0]
t = data[:,1]
u = data[:,2]

x_axis = x[0:J+1]
time = t[0:(J+1)*(N+1):J+1]
dt = time[1] - time[0]
Tmax = time[-1]

u.shape = (N+1,J+1)

files = []

for i in range(N+1):
    plt.cla()
    axes = plt.gca()
    axes.set_xlim([0,1])
    axes.set_ylim([0,110])
    plt.plot(x_axis, u[i,:], "k-")
    fname = '_tmp%03d.png' % i
    print('Saving frame', fname)
    plt.savefig(fname)
    files.append(fname)

print('Making movie Heat_equation.mpg - this may take a while')
fps = N / (Tmax * slowdown)
subprocess.call("mencoder 'mf://_tmp*.png' -mf type=png:fps=" + str(fps)
                + " -ovc lavc " +
                "-lavcopts vcodec=wmv2 -oac copy -o Heat_equation.mpg"
                , shell=True)
#TODO tidy up

#plt.show()

# cleanup
for fname in files:
    os.remove(fname)
