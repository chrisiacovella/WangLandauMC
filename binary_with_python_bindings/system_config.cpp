#include "WL.hpp"

configuration sys_config;

//Let us define all of the get and set functions we will need to perform the MC part of the code
extern void init_configuration(int N_particles_t)
{
    sys_config.init(N_particles_t);
}
//extern void set_configuration(std::vector < std::vector < double> > v_xyz_t, std::vector<int> v_type_t)
void set_configuration(coordlist_t v_xyz_t, std::vector<int> v_type_t)
{
    sys_config.v_xyz = v_xyz_t;
    sys_config.v_type = v_type_t;    
}

extern void get_configuration(std::vector < std::vector < double> >& v_xyz_t, std::vector<int>& v_type_t)
{
    v_xyz_t= sys_config.v_xyz;
    v_type_t = sys_config.v_type;
}

