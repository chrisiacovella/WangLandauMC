import WL

#define a few different containers using typedefs
xyz = WL.coordlist_t()
types = WL.types_t()
box = WL.coord_t()


#initialize the system
WL.init_system_binary(xyz, 110, 50, 0.5, box, types)

print "box lenght", box[0]
for i in range(0,10):
        print types[i], xyz[i][0], xyz[i][1], xyz[i][2]


print "\ncalculate some random gaussian numbers"
print WL.rand_gaussian()
print WL.rand_gaussian()
print WL.rand_gaussian()


print "calculate LJ potential at minimum"
#keep in mind this takes r2, not r
# not really functions we'll ever need to really call directly, but just for testing right now

r= 2.0**(1.0/6.0)
r2 = r*r
output =  WL.LJ_potential(r2, 0, 0)
print output

print "calculate LJ potential at overlapping condition, should return true"
r = 0.5
r2 = r*r
output =  WL.LJ_potential(r2, 0, 0)
print output

print "calculate the energy of the system"
pe = WL.calc_pe_global(xyz, types, box, 2.5)
print pe

print "calculate the energy of the system"
pe = WL.calc_pe_brute(xyz, types, box, 2.5)
print pe

print "calculate the energy of the system"
pe = WL.calc_pe(xyz, types, box, 2.5)
print pe

WL.print_xyz_file("test.xyz",  xyz, types, box)

mc = WL.metropolis_NVT()
mc.init_system(box, 12345, 0.25)
mc.init_translate(0.15,0.1, 0.25, 0.5)
mc.init_swap(50,0.1)
mc.init_nlist(2.5, 0.51)
mc.set_configuration(xyz, types)
mc.set_traj_filename("mc_traj_swap.xyz")
mc.set_thermo_filename("mc_thermo_swap.txt")
mc.setup()
mc.run(10000,1000, 1000, 100, True,False)
mc.run(50000,1000, 1000, 100, False,True)
