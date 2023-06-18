import WL
import math

xyz = WL.coordlist_t()
types = WL.types_t()
box = WL.coord_t()

number_density = 0.75
U_min = -200
U_max = 100

#initialize a binary system
WL.init_system_binary(xyz, 128, 64, 0.5, box, types)

#MC parameters
cutoff = 2.5
skin = 1.0
dx = skin/3.1
dx_max = skin/2.0001
dx_min = 0.001
dx_target = 0.5

#set up a simple metropolis MC run to converge the energy window
mNVT = WL.metropolis_NVT()
mNVT.init_system(box, 12345, 10.0)
mNVT.init_translate(dx, dx_min, dx_max, dx_target)
mNVT.set_configuration(xyz, types)
mNVT.init_nlist(cutoff,skin)
mNVT.set_traj_filename("mc_traj_swap.xyz")
mNVT.set_thermo_filename("mc_thermo_swap.txt")
mNVT.setup()

T_min = 0.1
T_max = 10.0

#converge the energy window
mNVT.converge_window(5000, 100, 10000, 5000, False, U_min, U_max, T_min, T_max)

#grab the configuration
mNVT.get_configuration(xyz, types)

print "system in range"

#WL parameters
U_min = -300
U_max = 100
U_bin = 2.0
flatness  = 0.7
N_min = 64
N_max = 124
N_bin = 4

#let's try out the WL part
mWL = WL.WL_2D()

#WL trial moves
#mWL.init_translate(dx, dx_min, dx_max, dx_target)
mWL.init_cswap(20,20,0.1)
mWL.init_swap(10, 10, 0.1)

#WL parameters
mWL.init_nlist(cutoff,skin)
mWL.init_system(box, 22345)
mWL.init_uniform_WL(U_min, U_max, U_bin, flatness, N_min, N_max, N_bin)
mWL.set_configuration(xyz, types)
mWL.set_traj_filename("WL_traj_swap.xyz")
mWL.set_thermo_filename("WL_thermo_swap.txt")
mWL.setup()

f = math.exp(1.0)

#mWL.calc_U_min(f, 10000, 1000, 1000, 1000, False, False)

isflat = False
count = 0
while f > math.exp(1e-4):
    while isflat == False:
        mNVT.set_temperature(10.0*WL.rand_double())
        mNVT.converge_window(5000, 100, 10000, 5000, False, U_min, U_max, T_min, T_max)
        print "current T: ", mNVT.T
        #mNVT.run(5000, 100, 1000, 1000, False, True)
        mNVT.get_configuration(xyz, types)
        mWL.set_configuration(xyz, types)
        
        if count == 10:
            isflat = mWL.run_step(f, False, True, True)
            count = 0
        else:
            isflat = mWL.run_step(f, False, True, False)
            count = count +1

        mWL.get_configuration(xyz, types)
        mNVT.set_configuration(xyz, types)

    mWL.clear_visited()
    isflat = False
    f = math.sqrt(f)
    print "******\n new f=",f
mWL.reset_configuration()


