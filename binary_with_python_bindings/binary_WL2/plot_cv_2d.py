import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl

import math 

#NP, U, V = np.loadtxt('visited_1.0000019074.txt').T
NP, U, V = np.loadtxt('ln_dos_temp.txt').T

max = len(NP)

print "max ", max

Umin = U.min()
Umax = U.max()
Vmin = V.min()
Vmax = V.max()

fig = plt.figure()
ax = fig.gca(projection='3d')
ax2 = fig.gca(projection='3d')

N_particles = 128


current_NP = NP[0]
next_NP = NP[0]
ii =0
while ii < max:
    UU = []
    SS = []
    NN = []
    T_min = 0.25
    T_max = 5.0
    T_delta = 0.01
    T_range = int((T_max-T_min)/T_delta)
    Z = []
    T = []
    E = []     #<E>
    E2 = []     #<E^2>
    E_2 = []    #<E>^2
    NA = []
    Cv = []
    current_T=T_min

    while next_NP == current_NP:
        SS.append(V[ii])
        UU.append(U[ii])
        NN.append(NP[ii])
        
        ii=ii+1
        if ii < max:
            next_NP = NP[ii]
        else:
            break;


    max2 = len(UU)
    for i in range(0, T_range):

        T.append(current_T)

        temp_Z=0
        temp_E = 0
        temp_E2 = 0
        tt = 0
        largest = -1000000000
        for j in range(0, max2):
            temp = SS[j]-UU[j]/current_T
            if temp > largest:
                largest = temp


        for j in range(0, max2):
            tt = math.exp((SS[j]-UU[j]/current_T)-largest)
            temp_Z += tt
            temp_E += tt*UU[j]
            temp_E2 += tt*UU[j]*UU[j]

        temp_E = temp_E/temp_Z;
        temp_E2 = temp_E2/(temp_Z);

        E.append(temp_E)
        E_2.append(temp_E*temp_E)
        E2.append(temp_E2)

        Cv_temp = (E2[i]-E_2[i])/(current_T*current_T*N_particles)
        Cv.append(Cv_temp)
        Z.append(temp_Z)
        NA.append(current_NP)
        #s1 = "%lf\t%lf\n" % (current_T, Cv_temp)
        #s2 = "%lf\t%lf\n" % (current_T, temp_E)
        #f_cv.write(str(s1))
        #f_tu.write(str(s2))
        current_T += T_delta


    ax.plot(NA, T, Cv)
#ax2.plot(NA, T, E)

    if ii < max:
        current_NP = NP[ii]
    else:
        break;

plt.show()


