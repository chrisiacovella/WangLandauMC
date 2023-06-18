/*
 *  energy.cpp
 *  simple_LJ
 *
 *  Created by Christopher Iacovella on 1/25/12.
 *  Copyright 2012. All rights reserved.
 *
 */

#include "WL.hpp"

//a simple wrapper that ruterns the potential energy total, using the other routine
double calc_pe_global(coordlist_t& x, types_t& types, coord_t L, double cutoff)
{
    double potential_energy_total;
    calc_pe_brute(x, types, L, &potential_energy_total, cutoff);
    return potential_energy_total;
}
//use this function to calculate system energy

void calc_pe_brute(coordlist_t& x, types_t& types, coord_t L, double *potential_energy_total, double cutoff)
{
	

	double cutoff2 = cutoff*cutoff;
		
	//precalculate 1/L
	double L_inv[3];
	for(int k=0; k<3; k++)
		L_inv[k] = 1.0/L[k];
		
	double potential_energy = 0.0;

	//brute force potential energy calculation, as it is done infrequently and does not need to be optimized
	for(int i=0; i<x.size(); i++)
	{
        
		for(int j=i+1; j<x.size(); j++)
		{
			
			//calculate the separation between two particles,
			//taking into account periodic boundary conditions
			double rij[3];
			double r2=0;
			for(int k=0; k<3; k++)
			{
				rij[k] =x[i][k] - x[j][k];
				double pbc  = L[k]*anint(rij[k]*L_inv[k]);
				rij[k] -= pbc;
				r2 += rij[k]*rij[k];
			}
            double pot_temp = 0;

            if(r2 < cutoff2)
            {
                bool overlap = LJ_potential(r2, types[i], types[j], &pot_temp);
                if(overlap)
                {
                    std::cerr << "warning particles " << i << "\t" << x[i][0] << "\t" << x[i][1] << "\t" << x[i][2] << "\tand\t" << j << "\t" << x[j][0] << "\t" << x[j][1] << "\t" << x[j][2] << " are overlapping in the configuration.\t" << sqrt(r2) << std::endl;
                }
                else
                    potential_energy += pot_temp;
			}
			
		}
		
	}
	*potential_energy_total = potential_energy;
	
}

bool calc_pe(coordlist_t& x, types_t& types, coord_t L, double *potential_energy_total, double cutoff)
{
	
    
	double cutoff2 = cutoff*cutoff;
    
	//precalculate 1/L
	double L_inv[3];
	for(int k=0; k<3; k++)
		L_inv[k] = 1.0/L[k];
    
	double potential_energy = 0.0;
    
	//brute force potential energy calculation, as it is done infrequently and does not need to be optimized
	for(int i=0; i<x.size(); i++)
	{
        
		for(int j=i+1; j<x.size(); j++)
		{
			
			//calculate the separation between two particles,
			//taking into account periodic boundary conditions
			double rij[3];
			double r2=0;
			for(int k=0; k<3; k++)
			{
				rij[k] =x[i][k] - x[j][k];
				double pbc  = L[k]*anint(rij[k]*L_inv[k]);
				rij[k] -= pbc;
				r2 += rij[k]*rij[k];
			}
            double pot_temp = 0;
            
            if(r2 < cutoff2)
            {
                bool overlap = LJ_potential(r2, types[i], types[j], &pot_temp);
                if(overlap)
                {
                    return true;
                }
                else
                    potential_energy += pot_temp;
			}
			
		}
		
	}
	*potential_energy_total = potential_energy;
	
    return false;
}
