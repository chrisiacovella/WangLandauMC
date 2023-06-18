//
//  neighbor.hpp
//  
//
//  Created by Chris Iacovella on 12/19/13.
//
//

#ifndef _neighbor_hpp
#define _neighbor_hpp
//neighborlist class holds a particle id and the initial position of a particle
class neighbor{
    
public:
	std::vector <int> member;
	coordlist_t x_old;
    coord_t dx;
	
};

void nsq_neighbor_init(coordlist_t& x, std::vector<neighbor> &nbr, double cutoff, double skin, coord_t L);
void nsq_neighbor_check_fast(coordlist_t& x, std::vector<neighbor> &nbr, double cutoff, double skin, coord_t L);
void nsq_neighbor_check(coordlist_t& x, std::vector<neighbor> &nbr, double cutoff, double skin, coord_t L);
void nsq_neighbor_rebuild(coordlist_t& x, std::vector<neighbor> &nbr, double cutoff, double skin, coord_t L);


#endif
