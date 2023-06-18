#ifndef _initialize_hpp
#define _initialize_hpp


extern void init_system(coordlist_t& x, int N, double density, coord_t& L, types_t& types);
extern void init_system_binary(coordlist_t & x, int N, int n_particles_0, double density, coord_t& L, types_t& types);

#endif
