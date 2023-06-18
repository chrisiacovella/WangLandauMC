
import WL
xyz = WL.coordlist_t()
xyz2 = WL.coordlist_t()
types = WL.types_t()
box = WL.coord_t()


for i in range (0,10):
    temp = WL.coord_t()
    temp.append(i)
    temp.append(1)
    temp.append(i*i)
    
    xyz.append(temp)



for i in range (0,10):
    print i, "\t", xyz[i][0], "\t", xyz[i][1], "\t", xyz[i][2]

for i in range (0,10):
    temp = WL.coord_t()
    temp.append(i)
    temp.append(i*i)
    temp.append(1)

    xyz[i] = temp


print "----"
for i in range (0,10):
    print i, "\t", xyz[i][0], "\t", xyz[i][1], "\t", xyz[i][2]

# or just run the "init routine, it will also set up box and types
WL.init_system_binary(xyz, 10, 5, 0.5, box, types)

for i in range (0,10):
    temp = WL.coord_t()
    temp.append(i)
    temp.append(i*i)
    temp.append(1)
    types[i] = 11
    xyz[i] = temp

print "----"
for i in range (0,10):
    print types[i], "\t", xyz[i][0], "\t", xyz[i][1], "\t", xyz[i][2]
