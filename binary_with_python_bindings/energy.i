%apply double *OUTPUT { double *potential_energy_total }
%{
#include "energy.hpp"
%}
%include "energy.hpp"


/*%{
double calc_pe_global(coordlist_t& x, types_t& types, coord_t L, double cutoff);
%}
double calc_pe_global(coordlist_t& x, types_t& types, coord_t L, double cutoff);
*/