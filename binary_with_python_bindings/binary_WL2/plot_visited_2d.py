import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl

import math 

#NP, U, V = np.loadtxt('ln_dos_1.0000152589.txt').T
#NP, U, V = np.loadtxt('ln_dos_1.0000000075.txt').T

#NP, U, V = np.loadtxt('ln_dos_1.2840254167.txt').T
#NP, U, V = np.loadtxt('visited_1.0317434075.txt').T
#NP, U, V = np.loadtxt('visited_2.7182818285.txt').T

#NP, U, V = np.loadtxt('visited_temp.txt').T

#NP, U, V = np.loadtxt('ln_dos_1.0000019074.txt').T
NP, U, V = np.loadtxt('visited_1.2840254167.txt').T

max = len(NP)

print "max ", max

Umin = U.min()
Umax = U.max()
Vmin = V.min()
Vmax = V.max()

fig = plt.figure()
ax = fig.gca(projection='3d')


current_NP = NP[0]
next_NP = NP[0]
i =0
while i < max:
    x = []
    y = []
    z = []
    while next_NP == current_NP:
        z.append(V[i])
        x.append(U[i])
        #        y.append(math.log(NP[i]))
        y.append(NP[i])
        i=i+1
        if i < max:
            next_NP = NP[i]
        else:
            break;

    ax.plot(x, y, z)
    if i < max:
        current_NP = NP[i]
    else:
        break;

plt.show()


