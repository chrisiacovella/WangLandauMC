/*
 *  metropolis.cpp
 *  simple_LJ
 *
 *  Created by Christopher Iacovella on 1/25/12.
 *  Copyright 2012. All rights reserved.
 *
 */

#include "WL.hpp"

metropolis_NVT::metropolis_NVT()
{
    translate_particles = false;
    system_initialized = false;
    nlist_initialized = false;
    configuration_initialized = false;
    swap_particles = false;
    setup_initialized = false;
    time_current = 0;
    xyzOut_initialized = false;
    thermoOut_initialized = false;
}
void metropolis_NVT::init_translate(double dx_t, double dx_min_t, double dx_max_t, double dx_target_prob_t)
{
    dx = dx_t;
    dx_max = dx_max_t;
    dx_min = dx_min_t;
    dx_target_prob = dx_target_prob_t;
    
    translate_particles = true;
}

void metropolis_NVT::init_system(coord_t L_t,  short unsigned int sseed_t, double T_t)
{
    L.clear();
    L.push_back(0); L.push_back(0); L.push_back(0);
    for(int k=0; k<3; k++)
        L[k] = L_t[k];
    
    sseed = sseed_t;
    srand48(sseed);
    T = T_t;
    beta = 1.0/T;
    
    system_initialized = true;
}

void metropolis_NVT::init_swap(int swap_t, double swap_target_prob_t)
{
    swap = swap_t;
    swap_target_prob = swap_target_prob_t;
    swap_particles = true;

    
}
void metropolis_NVT::init_nlist(double cutoff_t, double skin_t)
{
    
    cutoff = cutoff_t;
    skin = skin_t;
    
    nlist_initialized = true;
    
}
void metropolis_NVT::set_temperature(double T_t)
{
    T=T_t;
    beta = 1.0/T;
}

void metropolis_NVT::set_configuration(coordlist_t x_t, types_t types_array_t)
{
    N_particles = x_t.size();
    if(x.size() != x_t.size() )
    {
        x.clear();
        types.clear();
        
        for(int i=0; i< x_t.size(); i++)
        {
            x.push_back(x_t[i]);
            types.push_back(types_array_t[i]);
            
        }
    }
    else
    {
        for(int i=0; i< x_t.size(); i++)
        {
            x[i] = x_t[i];
            types[i] = types_array_t[i];
        }
    }
    configuration_initialized = true;
}

void metropolis_NVT::get_configuration(coordlist_t& x_t, types_t& types_array_t)
{
   if(x.size() != x_t.size())
   {
       x_t.clear();
       types_array_t.clear();
   
       for(int i=0; i< x.size(); i++)
       {
           x_t.push_back(x[i]);
           types_array_t.push_back(types[i]);
       }
   }
   else
   {
       for(int i=0; i< x_t.size(); i++)
       {
           x_t[i] = x[i];
           types_array_t[i] = types[i];
        }
    }
}
void metropolis_NVT::setup()
{
    assert(dx < skin/2.0);
    for(int k=0; k<3; k++)
        assert(cutoff < L[k]/2.0);
    assert(system_initialized);
    assert(nlist_initialized);
    assert(configuration_initialized);
    assert(x.size() > 0);
    assert(xyzOut_initialized);
    assert(thermoOut_initialized);
    nsq_neighbor_init(x, nbr, cutoff, skin, L);

    potential =0;
    calc_pe_brute(x, types, L, &potential,cutoff);
    setup_initialized = true;


}
void metropolis_NVT::set_traj_filename(const char* filename)
{
    strcpy(xyzFilename,(const char *) filename);
    xyzOut_initialized = true;
}
void metropolis_NVT::set_thermo_filename(const char* filename)
{
    strcpy(thermoFilename,(const char *) filename);
    thermoOut_initialized = true;
 
}

void metropolis_NVT::converge_window(int time_run, int time_U_update, int time_xyz_output, int time_thermo_output, bool adjust, double U_min, double U_max, double T_min, double T_max)
{
    
    
    run(time_run, time_U_update, time_xyz_output, time_thermo_output, adjust, false);
    

    bool out_of_range = true;
    while(out_of_range)
    {
        if(potential > U_max)
        {
            T /= 1.1;
            if(T < T_min)
                T = T_min;
        }
        else if (potential < U_min)
        {
            T *= 1.1;
            if(T > T_max)
            T = T_max;
        }
        else
        {
            out_of_range = false;
            std::cout << " converged." << std::endl;
            break;
        }
        std::cout<< "T\t" << T << std::endl;
        run(time_run, time_U_update, time_xyz_output, time_thermo_output, true, false);
        nsq_neighbor_check_fast(x, nbr, cutoff, skin, L);
        calc_pe_brute(x, types, L, &potential,cutoff);
        
     
    }
}


void metropolis_NVT::run(int time_run, int time_U_update, int time_xyz_output, int time_thermo_output, bool adjust, bool append_file)
{
    assert(setup_initialized == true);
    
    std::ofstream xyzOut;
    std::ofstream thermoOut;
    if(append_file)
    {
        xyzOut.open(xyzFilename,std::ios::app);
        thermoOut.open(thermoFilename,std::ios::app);
    }
    else
    {
        xyzOut.open(xyzFilename);
        thermoOut.open(thermoFilename);
    }
    
    beta = 1.0/T;
    
    //first calculate the system potential energy
    calc_pe_brute(x, types, L, &potential, cutoff);

    int n_accept_translate = 0;
    int n_accept_swap = 0;
    int since_last = 1;

    double translate_probability=0;
    double swap_probability=0;
    
    int time_adjust = 1000;
    if(time_thermo_output < time_adjust)
        time_adjust = time_thermo_output;
    
    for(int i=0; i<time_run; i++)
    {
        int n_accept_local;
        
        if(translate_particles)
        {
            n_accept_local = translate_NVT();
            n_accept_translate += n_accept_local;
        }
        //check to see if we need to rebuild the neighborlist
        nsq_neighbor_check_fast(x, nbr, cutoff, skin, L);
        if(swap_particles)
        {
            n_accept_local = swap_NVT(x, types, nbr, L, swap, &potential, beta);
            n_accept_swap += n_accept_local;
        }
        //various outputs
        if(time_current%time_U_update == 0)
        {
            calc_pe_brute(x, types, L, &potential, cutoff);

        }
        
        translate_probability = (double)n_accept_translate/(double)(since_last*N_particles);
        swap_probability = (double)n_accept_swap/(since_last*swap);
        
        if(time_current%time_thermo_output == 0)
        {


            thermoOut << time_current << "\t" <<  potential << "\t" <<  N_particles << "\t";
            if(translate_particles)
                    thermoOut <<  translate_probability << "\t";
            if(swap_particles)
                thermoOut <<  swap_probability << "\t";
            thermoOut << std::endl;

            std::cout << time_current << "\t" <<  potential << "\t" <<  N_particles << "\t";
            if(translate_particles)
                std::cout <<  translate_probability << "\t";
            if(swap_particles)
                std::cout <<  swap_probability << "\t";
            std::cout << std::endl;
            
            //print_thermo(thermoOut, time_current, rank, potential, N_particles, translate_probability, swap_probability, 0);
            //print_thermo(std::cout, time_current, rank, potential, N_particles, translate_probability, swap_probability, 0);
            
            
            if(translate_particles)
                std::cout << "\tdx:\t" << dx;
            if(swap_particles)
                std::cout << "\tn swaps:\t" << swap;
            std::cout<< std::endl;
            
        }
        if(time_current%time_adjust == 0 && adjust == true && time_current > 0)
        {

                  
            if(translate_probability < dx_target_prob  && time_current > time_thermo_output)
            {
                dx/=1.1;
                if(dx < dx_min)
                    dx = dx_min;
            }
            else if(translate_probability > dx_target_prob && time_current > time_thermo_output)
            {
                dx*=1.1;
                if(dx >dx_max)
                    dx = dx_max;
            }
            
            if(swap_probability < swap_target_prob && time_current > time_thermo_output)
                swap++;
            else if(swap_probability > swap_target_prob && time_current > time_thermo_output)
            {
                swap--;
                if(swap <1)
                    swap = 1;
            }
       
            n_accept_translate = 0;
            n_accept_swap = 0;
            since_last =0;
            
        }

        
        if(time_current%time_xyz_output == 0)
        {
            print_xyz(xyzOut, x, types, L);
        }
        
        time_current++;
        since_last++;
    }

}

int metropolis_NVT::swap_NVT(coordlist_t& x, std::vector<int>& types, std::vector<neighbor> &nbr, coord_t L, int n_swaps, double *potential_energy_total, double beta)
{
    beta = 1.0/T;
    double neg_beta = -beta;
    //precalculate 1/L
    
    //precalculate 1/L
    double L_inv[3];
    for(int k=0; k<3; k++)
        L_inv[k] = 1.0/L[k];
    
    
    double pot_temp = 0;
    
    int n_accept=0;
    std::vector<int> swapped;
    
    for(int i=0; i<x.size(); i++)
        swapped.push_back(0);
    
    for(int i=0; i<n_swaps; i++)
    {
        bool picking = true;
        
        int n1,t1;
        while(picking)
        {
            
            n1 = drand48()*x.size();
            t1 = types[n1];
            
            if(swapped[n1] < 1)
            {
                picking = false;
            }
        }
        picking = true;
        int n2, t2;
        while(picking)
        {
            n2 = drand48()*x.size();
            t2 = types[n2];
            
            
            if(swapped[n2] < 1)
            {
                
                if(t1 != t2)
                    picking = false;
            }
        }
        
        double potential_energy_new = 0.0;
        double potential_energy_old = 0.0;
        bool overlap1 = true;
        bool overlap2 = true;
        
        
        //calculate the potential energy of particle "i" with all of its neighbors "j"
        for(int j=1; j<nbr[n1].member.size(); j++)
        {
            int nbr_id = nbr[n1].member[j];
            
            //calculate the separation between two particles,
            //taking into account periodic boundary conditions
            double rij[3];
            double rij_new[3];
            double r2=0;
            
            for(int k=0; k<3; k++)
            {
                rij[k] = x[n1][k] - x[nbr_id][k];
                
                double pbc  = L[k]*anint(rij[k]*L_inv[k]);
                
                rij[k] -= pbc;
                r2 += rij[k]*rij[k];
            }
            
            pot_temp = 0;
            //assume no overlaps in this configuration
            LJ_potential(r2, t1, types[nbr_id], &pot_temp);
            potential_energy_old += pot_temp;
            
            int tt = types[nbr_id];
            if(nbr_id == n2)
                tt = t1;
            
            pot_temp = 0;
            overlap1 = LJ_potential(r2, t2, tt, &pot_temp);
            if(overlap1==true)
            {
                break;
            }
            potential_energy_new += pot_temp;
            
        }
        if(overlap1 == false)
        {
            for(int j=1; j<nbr[n2].member.size(); j++)
            {
                int nbr_id = nbr[n2].member[j];
                
                //calculate the separation between two particles,
                //taking into account periodic boundary conditions
                double rij[3];
                double rij_new[3];
                double r2=0;
                double r2_new=0;
                
                for(int k=0; k<3; k++)
                {
                    rij[k] = x[n2][k] - x[nbr_id][k];
                    
                    double pbc  = L[k]*anint(rij[k]*L_inv[k]);
                    
                    rij[k] -= pbc;
                    r2 += rij[k]*rij[k];
                }
                
                pot_temp = 0;
                //assume no overlaps in this configuration
                LJ_potential(r2, t2, types[nbr_id], &pot_temp);
                potential_energy_old += pot_temp;
                
                int tt = types[nbr_id];
                if(nbr_id == n1)
                    tt = t2;
                pot_temp = 0;

                overlap2 = LJ_potential(r2, t1, tt, &pot_temp);
                if(overlap2==true)
                {
                    break;
                }
                potential_energy_new += pot_temp;
                
            }
        }
        if(overlap1 == false && overlap2 == false)
        {
            //calculate the change in potential energy related to the displacement
            double delta_PE = potential_energy_new - potential_energy_old;
            
            
            //if we lower the potential energy always accept
            if(delta_PE < 0)
            {
                double pe_old = *potential_energy_total;
                double pe_new = *potential_energy_total + delta_PE;
                
                *potential_energy_total = pe_new;
                n_accept++;
                
                types[n1] = t2;
                types[n2] = t1;
//                coord_t temp_n1 = x[n1];
//                coord_t temp_n2 = x[n2];
//                neighbor nbr_n1 = nbr[n1];
//                neighbor nbr_n2= nbr[n2];
//                
//                x[n1] = temp_n2;
//                x[n2] = temp_n1;
//                nbr[n1] = nbr_n2;
//                nbr[n2] = nbr_n1;
                
                swapped[n1] = 1;
                swapped[n2] = 1;
            }
            else
            {
                double rand_value = drand48();
                if(exp(neg_beta*delta_PE) > rand_value)
                {
                    double pe_old = *potential_energy_total;
                    double pe_new = *potential_energy_total + delta_PE;
                    
                    *potential_energy_total = pe_new;
                    n_accept++;
                    
                    types[n1] = t2;
                    types[n2] = t1;
                    
                    swapped[n1] = 1;
                    swapped[n2] = 1;
                }
            }
            
        }
        
    }
    return n_accept;
    
}


int metropolis_NVT::translate_NVT()
{
    
    beta = 1.0/T;
    double neg_beta = -beta;
    //precalculate 1/L
    double L_inv[3];
    for(int k=0; k<3; k++)
        L_inv[k] = 1.0/L[k];
    
    
    double pot_temp = 0;
    
    int n_accept=0;
    for(int i=0; i<x.size(); i++)
    {
        
        double potential_energy_new = 0.0;
        double potential_energy_old = 0.0;
        bool overlap = false;
        
        //generate a random displacement of a particle
        //then check to see if we should accept this move
        coord_t new_x;
        coord_t new_dx;
        
        for( int k=0; k<3; k++)
        {
            double dx_temp = dx*(2.0*drand48()-1.0);
            new_x.push_back(x[i][k]+dx_temp);
            new_x[k] -= L[k]*anint(new_x[k]*L_inv[k]);
            
            new_dx.push_back(dx_temp);

        }
        //calculate the potential energy of particle "i" with all of its neighbors "j"
        for(int j=1; j<nbr[i].member.size(); j++)
        {
            int nbr_id = nbr[i].member[j];
            
            //calculate the separation between two particles,
            //taking into account periodic boundary conditions
            double rij[3];
            double rij_new[3];
            double r2=0;
            double r2_new=0;
            
            for(int k=0; k<3; k++)
            {
                rij[k] = x[i][k] - x[nbr_id][k];
                
                rij_new[k] = new_x[k] - x[nbr_id][k];
                
                double pbc  = L[k]*anint(rij[k]*L_inv[k]);
                double pbc_new  = L[k]*anint(rij_new[k]*L_inv[k]);
                
                rij[k] -= pbc;
                r2 += rij[k]*rij[k];
                
                rij_new[k] -= pbc_new;
                r2_new += rij_new[k]*rij_new[k];
            }
            
            pot_temp = 0;
            //assume no overlaps in this configuration
            LJ_potential(r2, types[i], types[nbr_id], &pot_temp);
            potential_energy_old += pot_temp;
            
            pot_temp = 0;
            overlap = LJ_potential(r2_new, types[i], types[nbr_id], &pot_temp);
            if(overlap==true)
            {
                break;
            }
            //std::cout << pot_temp << std::endl;
            potential_energy_new += pot_temp;
            
            
        }
        if(overlap == false)
        {
            //calculate the change in potential energy related to the displacement
            double delta_PE = potential_energy_new - potential_energy_old;
            
            //if we lower the potential energy always accept
            if(delta_PE <= 0)
            {
                potential += delta_PE;
                n_accept++;
                
                //assign the translated position to the main position array
                //making sure to apply PBC
                for( int k=0; k<3; k++)
                {
                    x[i][k] = new_x[k];
                    
                }
            }
            else
            {
                //if the potential energy increases, check to see if we should accept
                double rand_value = drand48();
                if(exp(neg_beta*delta_PE) > rand_value)
                {
                    potential += delta_PE;
                    n_accept++;
                    
                    //assign the translated position to the main position array
                    //making sure to apply PBC
                    for( int k=0; k<3; k++)
                    {
                        x[i][k] = new_x[k];
                        nbr[i].dx[k] += new_dx[k];
                        
                    }
                    
                }
            }
        }
        
    }
    return n_accept;
    
}

