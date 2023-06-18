import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math 

NP, U, S = np.loadtxt('S_all.txt').T

target_NP = 127

f_cv = open("Cv.txt", "w\n")
f_tu = open("TU.txt", "w\n")

N_particles = 128

max = len(U)
S_norm = S[0];
G = np.zeros(max)

    #for i in range(0,max):
#S[i] -= S_norm;

Umin = U.min()
Umax = U.max()
Smin = S.min()
Smax = S.max()


T_min = 0.25
T_max = 2.5
T_delta = 0.01
T_range = int((T_max-T_min)/T_delta)
Z = []
T = []
E = []     #<E>
E2 = []     #<E^2>
E_2 = []    #<E>^2

Cv = []
current_T=T_min
for i in range(0, T_range):
    
    T.append(current_T)
    
    temp_Z=0
    temp_E = 0
    temp_E2 = 0
    tt = 0
    largest = -1000000000
    for j in range(0, max):
        temp = S[j]-U[j]/current_T
        if temp > largest:
            largest = temp


    for j in range(0, max):

        if NP[j] == target_NP:
            tt = math.exp((S[j]-U[j]/current_T)-largest)
            temp_Z += tt
            temp_E += tt*U[j]
            temp_E2 += tt*U[j]*U[j]
    
    temp_E = temp_E/temp_Z;
    temp_E2 = temp_E2/(temp_Z);

    E.append(temp_E)
    E_2.append(temp_E*temp_E)
    E2.append(temp_E2)

    Cv_temp = (E2[i]-E_2[i])/(current_T*current_T*N_particles)
    Cv.append(Cv_temp)
    Z.append(temp_Z)
    s1 = "%lf\t%lf\n" % (current_T, Cv_temp)
    s2 = "%lf\t%lf\n" % (current_T, temp_E)
    f_cv.write(str(s1))
    f_tu.write(str(s2))
    current_T += T_delta


#plt.plot(U, S)

f_cv.close()
f_tu.close()
plt.plot(T,Cv)
#plt.axis([0.5,2,0,10])
plt.savefig('Cv_T.png')

plt.show()

plt.plot(E,T)
plt.savefig('T_E.png')
plt.show()

