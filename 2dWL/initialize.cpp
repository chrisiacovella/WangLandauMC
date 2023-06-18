/*
 *  initialize.cpp
 *  simple_LJ
 *
 *  Created by Christopher Iacovella on 8/2/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "main.h"


//randomly place N particles in the box at a specific density
void init_system(coordlist_t& x, int N, double density, double *L, std::vector<int>& types)
{
	int i = 0;
	double cut2 = 1.0*1.0;
	
	L[0] = pow(N/density, (1.0/3.0));
	L[1] = L[0];
	L[2] = L[0];
	double box = L[0]; //since cubic
	
    //precalculate 1/L
	double L_inv[3];
	for(int k=0; k<3; k++)
		L_inv[k] = 1.0/L[k];
    
	
	while(i < N)
	{
        coord_t xi;
        xi.clear();
        for(int k=0; k<3; k++)
            xi.push_back(drand48() * L[k] - 0.5 * L[k]);
		
		bool overlap = false;
		for(int j=0; j<i; j++)
		{
            
            double rij[3];
			double r2=0;
			for(int k=0; k<3; k++)
			{
				rij[k] =xi[k] - x[j][k];
				double pbc  = L[k]*anint(rij[k]*L_inv[k]);
				rij[k] -= pbc;
				r2 += rij[k]*rij[k];
			}
            
			if(r2 < cut2)
			{
				overlap = true;
				break;
			}
		}
		
		if(overlap == false)
		{
			x.push_back(xi);
            types.push_back(0);
			i++;
            
		}
	}
}

//randomly place N particles in the box at a specific density
void init_system_binary(coordlist_t& x, int N, int n_particles_0, double density, double *L, std::vector<int>& types)
{
	int i = 0;
	double cut2 = 0.88*0.88;
	
	L[0] = pow(N/density, (1.0/3.0));
	L[1] = L[0];
	L[2] = L[0];
	double box = L[0]; //since cubic
	
    //precalculate 1/L
	double L_inv[3];
	for(int k=0; k<3; k++)
		L_inv[k] = 1.0/L[k];
    
    int N0 = n_particles_0;
	
	while(i < N)
	{
        coord_t xi;
        xi.clear();
        for(int k=0; k<3; k++)
            xi.push_back(drand48() * L[k] - 0.5 * L[k]);
		
		bool overlap = false;
		for(int j=0; j<i; j++)
		{
            
            double rij[3];
			double r2=0;
			for(int k=0; k<3; k++)
			{
				rij[k] =xi[k] - x[j][k];
				double pbc  = L[k]*anint(rij[k]*L_inv[k]);
				rij[k] -= pbc;
				r2 += rij[k]*rij[k];
			}
            
			if(r2 < cut2)
			{
				overlap = true;
				break;
			}
		}
		
		if(overlap == false)
		{
			x.push_back(xi);
            if(i<N0)
                types.push_back(0);
            else
                types.push_back(1);

			i++;
            
		}
	}
}

//generates a pseudo random gaussian distribution
double rand_gaussian()
{
	double first,v1,v2,rsq,fac;
	static int save = 0;
	static double second;
	
	if (!save) {
		int again = 1;
		while (again) {
			v1 = 2.0*drand48()-1.0;
			v2 = 2.0*drand48()-1.0;
			rsq = v1*v1 + v2*v2;
			if (rsq < 1.0 && rsq != 0.0) again = 0;
		}
		fac = sqrt(-2.0*log(rsq)/rsq);
		second = v1*fac;
		first = v2*fac;
		save = 1;
	} 
	else {
		first = second;
		save = 0;
	}
	return first;
}

void check_cutoff(double cutoff, double skin, double *L)
{
    
    //make sure we didn't do anything bad
    if(cutoff >= L[0]/2.0)
        std::cerr << "cutoff can't be larger than half the box. cutoff: "<< cutoff <<"\t L[0]/2.0\t" << L[0]/2.0 << std::endl;
    if(cutoff >= L[1]/2.0)
        std::cerr << "cutoff can't be larger than half the box. cutoff: "<<cutoff <<"\t L[1]/2.0\t" << L[1]/2.0 << std::endl;
    if(cutoff >= L[2]/2.0)
        std::cerr << "cutoff can't be larger than half the box. cutoff: "<<cutoff <<"\t L[2]/2.0\t" << L[2]/2.0 << std::endl;
    
    assert(cutoff < L[0]/2.0);
    assert(cutoff < L[1]/2.0);
    assert(cutoff < L[2]/2.0);
    
    if(cutoff+skin >= L[0]/2.0)
        std::cerr << "cutoff+skin can't be larger than half the box. cutoff: "<<cutoff <<"\t skin " << skin << "\t L[0]/2\t" << L[0]/2.0 << std::endl;
    if(cutoff+skin >= L[1]/2.0)
        std::cerr << "cutoff+skin can't be larger than half the box. cutoff: "<<cutoff <<"\t skin " << skin << "\t L[1]/2\t" << L[1]/2.0 << std::endl;
    if(cutoff+skin >= L[2]/2.0)
        std::cerr << "cutoff+skin can't be larger than half the box. cutoff: "<<cutoff <<"\t skin " << skin << "\t L[2]/2\t" << L[2]/2.0 << std::endl;
    
}



