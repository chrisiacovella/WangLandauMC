/*
 *  energy.cpp
 *  simple_LJ
 *
 *  Created by Christopher Iacovella on 1/25/12.
 *  Copyright 2012. All rights reserved.
 *
 */

#include "main.h"


//let's just define the interaction energy in one spot to make life easier
//will make life easier when we have more than one type of particle
//inline bool LJ_potential(double r2, int p1_type, int p2_type, double *pe)
//{
//    int pair_combo = p1_type*10+p2_type*10;
//        
//    //precalculate a bunch of terms for efficiency
//    double r_inv_p2 = 1.0/r2;                       // 1/r^2
//    double r_inv_p6 = r_inv_p2*r_inv_p2*r_inv_p2;	// 1/r^6
//    double r_inv_p12 = r_inv_p6*r_inv_p6;           // 1/r^12
//    
//    double pot = 0;
//    
//    
//    //6.25 = 2.5*2.5
//    //0.7225 = 0.85*0.85
//    if(pair_combo == 0)
//    {
//        if(r2 > 0.7225)
//        {
//            if(r2 < 6.25)
//            {
//                double sigma_p12_4eps = 4.0; //4.0*epsilon*sigma^12
//                double sigma_p6_4eps = 4.0; //4.0*epsilon*sigma^6
//                
//                pot = sigma_p12_4eps*r_inv_p12 - sigma_p6_4eps*r_inv_p6;
//                *pe = pot;
//                return false;
//
//            }
//
//        }
//        else             //overlapping
//        {
//            *pe = pot;
//            return true;
//        }
//    }
//    else if(pair_combo == 10)
//    {
//        if(r2 > 0.7225)
//        {
//            if(r2 < 6.25)
//            {
//                double sigma_p12_4eps = 2.0; //4.0*epsilon*sigma^12, epsilon = 0.5
//                double sigma_p6_4eps = 2.0; //4.0*epsilon*sigma^6, epsilon = 0.5
//                
//                pot = sigma_p12_4eps*r_inv_p12 - sigma_p6_4eps*r_inv_p6;
//                *pe = pot;
//                return false;
//                
//            }
//            
//        }
//        else             //overlapping
//        {
//            *pe = pot;
//            return true;
//        }
//    }
//    else if(pair_combo == 20)
//    {
//        if(r2 > 0.7225)
//        {
//            if(r2 < 6.25)
//            {
//                double sigma_p12_4eps = 4.0; //4.0*epsilon*sigma^12, epsilon = 1.0
//                double sigma_p6_4eps = 4.0; //4.0*epsilon*sigma^6, epsilon = 1.0
//                
//                pot = sigma_p12_4eps*r_inv_p12 - sigma_p6_4eps*r_inv_p6;
//                *pe = pot;
//                return false;
//                
//            }
//            
//        }
//        else             //overlapping
//        {
//            *pe = pot;
//            return true;
//        }
//    }
//
//    return false;
//}







//use this function to calculate system energy
void calc_pe_brute(coordlist_t& x, std::vector<int>& types, double *L, double *potential_energy_total, double cutoff)
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
                    std::cerr << "warning particles " << i << "\t" << x[i] << "\tand\t" << j << "\t" << x[j] << " are overlapping in the configuration.\t" << sqrt(r2) << std::endl;
                }
                else
                    potential_energy += pot_temp;
			}
			
		}
		
	}
	*potential_energy_total = potential_energy;
	
}
