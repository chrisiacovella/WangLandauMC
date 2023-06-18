//
//  energy.hpp
//  
//
//  Created by Chris Iacovella on 12/19/13.
//
//

#ifndef _energy_hpp
#define _energy_hpp
double calc_pe_global(coordlist_t& x, types_t& types, coord_t L, double cutoff);
void calc_pe_brute(coordlist_t& x, types_t& types, coord_t L, double *potential_energy_total, double cutoff);
bool calc_pe(coordlist_t& x, types_t& types, coord_t L, double *potential_energy_total, double cutoff);

#endif
