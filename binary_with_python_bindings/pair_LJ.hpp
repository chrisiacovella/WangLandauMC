/*
 *  pair_LJ.hpp
 *  simple_LJ
 *
 *  Created by Christopher Iacovella on 1/25/12.
 *  Copyright 2012. All rights reserved.
 *
 */

#ifndef _pair_LJ_hpp
#define _pair_LJ_hpp

//let's just define the interaction energy in one spot to make life easier
//will make life easier when we have more than one type of particle
//we will also define the swap function, so that you can change a particle type
//in a composition change
inline int swap_type(int p1_type)
{
    if(p1_type == 0)
    {
        return 1;
    }
    else
        return 0;
}



inline bool LJ_potential(double r2, int p1_type, int p2_type, double *pe)
{
    int pair_combo = p1_type*10+p2_type*10;
        
    //precalculate a bunch of terms for efficiency
    double r_inv_p2 = 1.0/r2;                       // 1/r^2
    double r_inv_p6 = r_inv_p2*r_inv_p2*r_inv_p2;	// 1/r^6
    double r_inv_p12 = r_inv_p6*r_inv_p6;           // 1/r^12
    
    double pot = 0;
    
    
    //6.25 = 2.5*2.5
    //0.7225 = 0.85*0.85
    if(pair_combo == 0)
    {
        if(r2 > 0.7225)
        {
            if(r2 < 6.25)
            {
                double sigma_p12_4eps = 4.0; //4.0*epsilon*sigma^12
                double sigma_p6_4eps = 4.0; //4.0*epsilon*sigma^6
                
                pot = sigma_p12_4eps*r_inv_p12 - sigma_p6_4eps*r_inv_p6;
                *pe = pot;
                return false;

            }

        }
        else             //overlapping
        {
            *pe = pot;
            return true;
        }
    }
    else if(pair_combo == 10)
    {
        if(r2 > 0.7225)
        {
            if(r2 < 6.25)
            {
                double sigma_p12_4eps = 0.4; //4.0*epsilon*sigma^12, epsilon = 0.1
                double sigma_p6_4eps = 0.4; //4.0*epsilon*sigma^6, epsilon = 0.1
                
                pot = sigma_p12_4eps*r_inv_p12 - sigma_p6_4eps*r_inv_p6;
                *pe = pot;
                return false;
                
            }
            
        }
        else             //overlapping
        {
            *pe = pot;
            return true;
        }
    }
    else if(pair_combo == 20)
    {
        if(r2 > 0.7225)
        {
            if(r2 < 6.25)
            {
                double sigma_p12_4eps = 4.0; //4.0*epsilon*sigma^12, epsilon = 1.0
                double sigma_p6_4eps = 4.0; //4.0*epsilon*sigma^6, epsilon = 1.0
                
                pot = sigma_p12_4eps*r_inv_p12 - sigma_p6_4eps*r_inv_p6;
                *pe = pot;
                return false;
                
            }
            
        }
        else             //overlapping
        {
            *pe = pot;
            return true;
        }
    }

    return false;
}
#endif