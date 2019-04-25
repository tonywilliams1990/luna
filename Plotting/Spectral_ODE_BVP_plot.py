import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
import argparse

lines = ["-","--","-.",":"]
linecycler = cycle(lines)

CLI=argparse.ArgumentParser()
CLI.add_argument(
  "--n",  # name on the CLI - drop the `--` for positional/required parameters
  nargs="*",  # 0 or more values expected => creates a list
  type=int,
  default=[],  # default if nothing is provided
)

# parse the command line
args = CLI.parse_args()
ns = args.n

plt.figure()

for i in range(len(ns)):
    data = np.loadtxt( "./DATA/Spectral_ODE_BVP.dat" )
    x = data[:,0]
    u_exact = data[:,1]
    u_spectral = data[:,i+2]

    plt.plot( x, u_exact, 'r', label="exact" )
    line = next(linecycler)
    plt.plot( x, u_spectral, c='k', linestyle=line, label="n = " + str( ns[ i ] ) )

plt.legend()
plt.show()
