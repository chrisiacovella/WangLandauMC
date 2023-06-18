/*
 *  wanglandau.hpp
 *  simple_LJ
 *
 *  Created by Christopher Iacovella
 *  Copyright 2012. All rights reserved.
 *
 */


//This is a class to handle standard metropolis MC
//This class will contain all the basic routines to minimize the amount of "stuff" in the main function

#include "WL.hpp"

WL_2D::WL_2D()
{
    translate_particles = false;
    system_initialized = false;
    nlist_initialized = false;
    configuration_initialized = false;
    swap_particles = false;
    cswap_particles = false;
    WL_initialized = false;
    setup_initialized = false;
    vol_change = false;
    xyzOut_initialized =false;
    thermoOut_initialized =false;
    loggerOut_initialized = false;
    nn=0;
    time_current = 0;
}

void WL_2D::set_traj_filename(const char* filename)
{
    strcpy(xyzFilename,(const char *) filename);
    xyzOut_initialized = true;
}
void WL_2D::set_thermo_filename(const char* filename)
{
    strcpy(thermoFilename,(const char *) filename);
    thermoOut_initialized = true;
    
}

void WL_2D::init_uniform_WL(double U_min_t, double U_max_t, double U_bin_t, double flatness_t, double M_min_t, double M_max_t, double M_bin_t)
{

    flatness = flatness_t;
    
    M_min = M_min_t;
    M_max = M_max_t;
    M_bin = M_bin_t;
    
    U_max_global = U_max_t;
    U_bin_global = U_bin_t;
    visited.init_uniform(U_min_t, U_max_t, U_bin_t, M_min, M_max, M_bin);
    
    ln_dos.init_uniform(U_min_t, U_max_t, U_bin_t, M_min, M_max, M_bin);
    visited.check();
    ln_dos.check();
    WL_initialized = true;
    
}

void WL_2D::init_from_params_WL(const char* filename, double flatness_t)
{
   
    flatness = flatness_t;
    
    visited.read_histogram_params(filename);
    
    ln_dos.read_histogram_params(filename);
    visited.check();
    ln_dos.check();
    WL_initialized = true;
    
}

void WL_2D::init_ln_dos_from_file(const char* filename_ln_dos)
{
    
    if(WL_initialized == true)
    {
        ln_dos.read_from_file(filename_ln_dos, 1);
    }
    else
    {
        std::cerr<< "please define energy range and binsize first in init_WL" << std::endl;
        assert(WL_initialized==true);
    }
}

void WL_2D::init_y(double M_min_t, double M_max_t, double M_bin_t)
{
    M_min = M_min_t;
    M_max = M_max_t;
    M_bin = M_bin_t;
    visited.init_y(M_min, M_max, M_bin);
    ln_dos.init_y(M_min, M_max, M_bin);
    
}

void WL_2D::push_back_x(double U_min_t, double U_max_t, double U_bin_t)
{
    visited.push_back_x(U_min_t, U_max_t, U_bin_t);
    ln_dos.push_back_x(U_min_t, U_max_t, U_bin_t);
}

void WL_2D::setup_histograms(double flatness_t)
{
    visited.setup();
    ln_dos.setup();
    visited.check();
    ln_dos.check();
    
    flatness = flatness_t;
    WL_initialized = true;

}

void WL_2D::init_translate(double dx_t, double dx_min_t, double dx_max_t, double dx_target_prob_t)
{
    dx = dx_t;
    dx_max = dx_max_t;
    dx_min = dx_min_t;
    dx_target_prob = dx_target_prob_t;
    
    translate_particles = true;
}

void WL_2D::init_system(coord_t L_t,  short unsigned int sseed_t)
{
    L.clear();
    L_init.clear();
    for(int k=0; k<3; k++)
    {
        L.push_back(L_t[k]);
        L_init.push_back(L_t[k]);
    }
    sseed = sseed_t;
    system_initialized = true;
}

void WL_2D::init_swap(int swap_t, int swap_max_t, double swap_target_prob_t)
{
    swap = swap_t;
    swap_max = swap_max_t;
    swap_target_prob = swap_target_prob_t;
    swap_particles = true;
    
    
}
void WL_2D::init_cswap(int cswap_t, int cswap_max_t, double cswap_target_prob_t)
{
    if(vol_change == false)
    {
        cswap = cswap_t;
        cswap_max = cswap_max_t;
        cswap_target_prob= cswap_target_prob_t;
        cswap_particles = true;
    }
    else
    {
        std::cerr << "the 2nd WL dimension can either be a volume change or a composition change, not both at this time" << std::endl;
        assert(vol_change == false);
    }
    
}

void WL_2D::init_vol_change(int vol_t, int vol_max_t, double vol_target_prob_t, double ln_vol_change_max_t)
{
    if(cswap_particles == false)
    {
        vol = vol_t;
        vol_max = vol_max_t;
        vol_target_prob = vol_target_prob_t;
        ln_vol_change_max = ln_vol_change_max_t;
        vol_change = true;
    }
    else
    {
        std::cerr << "the 2nd WL dimension can either be a volume change or a composition change, not both at this time" << std::endl;
        assert(cswap_particles == false);
    }
    
}


void WL_2D::init_nlist(double cutoff_t, double skin_t)
{
    cutoff = cutoff_t;
    skin = skin_t;
    
    nlist_initialized = true;
    
}
void WL_2D::set_configuration(coordlist_t x_t, types_t types_array_t)
{
    N_particles = x_t.size();
    if(x.size() != x_t.size() )
    {
        x.clear();
        x_init.clear();
        
        types.clear();
        types_init.clear();
        
        for(int i=0; i< x_t.size(); i++)
        {
            x.push_back(x_t[i]);
            x_init.push_back(x_t[i]);
            
            types.push_back(types_array_t[i]);
            types_init.push_back(types_array_t[i]);
            
        }
    }
    else
    {
        for(int i=0; i< x_t.size(); i++)
        {
            x[i] = x_t[i];
            x_init[i] = x_t[i];
            types[i] = types_array_t[i];
            types_init[i] = types_array_t[i];
        }
    }
    configuration_initialized = true;
}
void WL_2D::reset_configuration()
{
    assert(x.size() == x_init.size() );
    
    for(int i=0; i<x.size(); i++)
    {
        x[i] = x_init[i];
        types[i] = types_init[i];
    }
    for(int k=0; k<3; k++)
        L[k] = L_init[k];
    
    nsq_neighbor_rebuild(x, nbr, cutoff, skin, L);
    
    if(cswap_particles == true)
        count_type0();
    
}



void WL_2D::get_configuration(coordlist_t& x_t, types_t& types_array_t)
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
void WL_2D::setup()
{
    assert(dx < skin/2.0);
    for(int k=0; k<3; k++)
        assert(cutoff <= L[k]/2.0);
    assert(system_initialized);
    assert(nlist_initialized);
    assert(configuration_initialized);
    assert(WL_initialized);
    assert(x.size() > 0);
    assert(xyzOut_initialized);
    assert(thermoOut_initialized);
    nsq_neighbor_init(x, nbr, cutoff, skin, L);
    
    potential =0;
    calc_pe_brute(x, types, L, &potential, cutoff);
    if(translate_particles)
        nn+=x.size();
    

    sprintf(loggerFilename, "log.txt");
    setup_initialized = true;
    
    n_accept_translate = 0;
    n_accept_swap = 0;
    n_accept_cswap = 0;
    n_accept_vol = 0;
    iterations = 1;
    
    
}

void WL_2D::calc_U_min(double f_t, int run_time, int time_U_update, int time_check_flat, int time_thermo_output, bool adjust,bool append_file)
{
    
    time_current = 0;
    since_last=0;
    assert(setup_initialized == true);
    
    std::ofstream loggerOut;
    
    
    if(append_file)
    {
        loggerOut.open(loggerFilename,std::ios::app);
    }
    else
    {
        loggerOut.open(loggerFilename);
    }

    
    f = f_t;
    ln_f = log(f);
    
    ln_dos_old = ln_dos;
    //first calculate the system potential energy
    potential =0;
    calc_pe_brute(x, types, L, &potential, cutoff);
    
    n_accept_translate = 0;
    n_accept_swap = 0;
    n_accept_cswap = 0;
    n_accept_vol = 0;
    
    int time_since=0;
    //make sure we clear out the visited histogram
    visited.clear_all();
    
    bool is_flat = false;
    
    //we need to set up M_current
    if(cswap_particles == true)
        count_type0();
    
    if(vol_change == true)
        calc_vol();
    
    int time_local = 0;
    while(time_current < run_time)
    {
        if(is_flat == true)
            break;
        int n_accept_local;
        
        
        if(translate_particles)
        {
            n_accept_local = translate_WL();
            n_accept_translate += n_accept_local;
        }
        
        nsq_neighbor_check_fast(x, nbr, cutoff, skin, L);
        
        //check to see if we need to rebuild the neighborlist
        if(swap_particles)
        {
            n_accept_local = swap_WL();
            n_accept_swap += n_accept_local;
        }
        
        if(cswap_particles)
        {
            n_accept_local = cswap_WL();
            n_accept_cswap += n_accept_local;
        }
        
        if(vol_change)
        {
            n_accept_local = vol_WL();
            n_accept_vol += n_accept_local;
        }
        
        //various outputs
        if(time_current%time_U_update == 0)
        {
            potential =0;
            calc_pe_brute(x, types, L, &potential, cutoff);
        }
        double translate_probability = (double)n_accept_translate/(double)(since_last*x.size());
        cswap = x.size();
        double cswap_probability = (double)n_accept_cswap/(since_last*cswap);
        double swap_probability = (double)n_accept_swap/(since_last*swap);
        double vol_probability = (double)n_accept_vol/(since_last*vol);
        
        if(time_current%time_thermo_output == 0)
        {
            std::cout << "M_current: " << M_current << std::endl;
            

            
            std::cout << time_current << "\t" << potential << "\t";
            if(translate_particles)
            {
                std::cout << "trans:\t"<< translate_probability << "\t";
            }
            if(swap_particles)
            {
                std::cout << "swap:\t" << swap_probability << "\t";
            }
            if(cswap_particles)
            {
                std::cout << "cswap:\t" << cswap_probability << "\t";
            }
            if(vol_change)
            {
                std::cout << "vol:\t" << vol_probability << "\t";
            }
            std::cout << std::endl;
            
      
            //print_thermo(std::cout, time_current, rank, potential, N_particles, translate_probability, cswap_probability, vol_probability);
            
            if(translate_particles)
                std::cout << "\tdx:\t" << dx;
            if(cswap_particles)
                std::cout << "\tn swaps:\t" << cswap;
            std::cout<< std::endl;
            
            
            loggerOut << time_current << "\t" << potential << std::endl;
            
            char f1[1000];
            sprintf(f1, "visited_temp_converge.txt");
            visited.print_histogram(f1);
            
            
            if(adjust == true)
            {
                if(translate_particles)
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
                    n_accept_translate = 0;
                    translate_probability  = 0;
                    
                }
                if(cswap_particles)
                {
                    if(cswap_probability < cswap_target_prob && time_current > time_thermo_output)
                    {
                        cswap++;
                        if(cswap > cswap_max)
                            cswap = cswap_max;
                    }
                    else if(cswap_probability > cswap_target_prob && time_current > time_thermo_output)
                    {
                        cswap--;
                        if(cswap < 1)
                        {
                            cswap = 1;
                        }
                    }
                    n_accept_cswap = 0;
                    cswap_probability = 0;
                    
                }
                n_accept_translate = 0;
                n_accept_cswap = 0;
                n_accept_swap = 0;
                n_accept_vol = 0;
                since_last=0;
            }
            
        }
        
        
        if(time_current%time_check_flat == 0 && time_current > time_check_flat)
        {
            
            double avg_flat=0, min_flat=0, max_flat=0, metric=0;
            bool check_flat =  check_flatness(visited, flatness, &avg_flat, &min_flat, &max_flat, &metric);

                
            std::cout << std::setprecision(5) << "flatness min\t" << min_flat << "\tmax\t" << max_flat << "\tm\t" << metric << std::endl;
            if(check_flat == true)
            {
                std::cout << "histogram converged with " << flatness << " for f = " << std::setprecision(32) << f << std::endl;
                
            }
            
            if(check_flat == true)
            {
                is_flat = true;
            }
            
        }
        time_since++;
        time_current++;
        time_local++;
        since_last++;
        
    }
    
    std::vector<double> U_min_array;
    
    U_min_array = visited.calculate_min();
    visited.reset();
    ln_dos.reset();
    
    
    visited.init_y(M_min, M_max, M_bin);
    ln_dos.init_y(M_min, M_max, M_bin);
    
    
    for(int i=0; i<U_min_array.size(); i++)
    {
        visited.push_back_x(U_min_array[i], U_max_global, U_bin_global);
        ln_dos.push_back_x(U_min_array[i], U_max_global, U_bin_global);
        std::cout << "bin\t" << visited.container[i].M << "\tU min:\t" << U_min_array[i] << std::endl;

    }
    visited.setup();
    ln_dos.setup();

    visited.print_params("histogram_params.txt");
    
}
void WL_2D::clear_visited()
{
    visited.clear_all();
}
void WL_2D::run(double f_t, int time_U_update, int time_check_flat, int time_xyz_output, int time_thermo_output, bool adjust, bool append_file)
{
    
    time_current = 0;
    since_last = 0;
    assert(setup_initialized == true);
    
    std::ofstream xyzOut;
    std::ofstream thermoOut;
    std::ofstream loggerOut;
    
    
    if(append_file)
    {
        xyzOut.open(xyzFilename,std::ios::app);
        thermoOut.open(thermoFilename,std::ios::app);
        loggerOut.open(loggerFilename,std::ios::app);
    }
    else
    {
        xyzOut.open(xyzFilename);
        thermoOut.open(thermoFilename);
        loggerOut.open(loggerFilename);

    }
    
    
    f = f_t;
    ln_f = log(f);
    
    ln_dos_old = ln_dos;
    //first calculate the system potential energy
    potential =0;
    calc_pe_brute(x, types, L, &potential, cutoff);
    
    int n_accept_translate = 0;
    int n_accept_swap = 0;
    int n_accept_cswap = 0;
    int n_accept_vol = 0;
    cswap = x.size();

    int time_since=0;
    //make sure we clear out the visited histogram
    visited.clear_all();

    bool is_flat = false;
    
    //we need to set up M_current
    if(cswap_particles == true)
        count_type0();
    
    if(vol_change == true)
        calc_vol();

    int time_local = 0;
    while(is_flat == false)
    {
        

        int n_accept_local;

        
        if(translate_particles)
        {
            n_accept_local = translate_WL();
            n_accept_translate += n_accept_local;
        }

        nsq_neighbor_check_fast(x, nbr, cutoff, skin, L);
        
        //check to see if we need to rebuild the neighborlist
        if(swap_particles)
        {
            n_accept_local = swap_WL();
            n_accept_swap += n_accept_local;
        }

        if(cswap_particles)
        {
            n_accept_local = cswap_WL();
            n_accept_cswap += n_accept_local;
        }

        if(vol_change)
        {
            n_accept_local = vol_WL();
            n_accept_vol += n_accept_local;
            
            nsq_neighbor_rebuild(x, nbr, cutoff, skin, L);
        }
        
        //various outputs
        if(time_current%time_U_update == 0)
        {
            potential =0;
            calc_pe_brute(x, types, L, &potential, cutoff);
        }
        double translate_probability = (double)n_accept_translate/(double)(since_last*x.size());
        double cswap_probability = (double)n_accept_cswap/(since_last*cswap);
        double swap_probability = (double)n_accept_swap/(since_last*swap);
        double vol_probability = (double)n_accept_vol/(since_last*vol);
        
        if(time_current%time_thermo_output == 0)
        {
            std::cout << "M_current: " << M_current << std::endl;

            
            thermoOut << time_current << "\t" << potential << "\t";
            std::cout << time_current << "\t" << potential << "\t";
            if(translate_particles)
            {
                thermoOut << translate_probability << "\t";
                std::cout << "trans:\t"<< translate_probability << "\t";
            }
            if(swap_particles)
            {
                thermoOut << swap_probability << "\t";
                std::cout << "swap:\t" << swap_probability << "\t";
            }
            if(cswap_particles)
            {
                thermoOut << cswap_probability << "\t";
                std::cout << "cswap:\t" << cswap_probability << "\t";
            }
            if(vol_change)
            {
                thermoOut << vol_probability << "\t";
                std::cout << "vol:\t" << vol_probability << "\t";
            }
            thermoOut << std::endl;
            std::cout << std::endl;

            
            //print_thermo(thermoOut, time_current, rank, potential, N_particles, translate_probability, cswap_probability, vol_probability);
            //print_thermo(std::cout, time_current, rank, potential, N_particles, translate_probability, cswap_probability, vol_probability);
            
            
            if(translate_particles)
                std::cout << "\tdx:\t" << dx;
            if(cswap_particles)
                std::cout << "\tn swaps:\t" << cswap;
            std::cout<< std::endl;
            

            loggerOut << time_current << "\t" << potential << std::endl;
            
            char f1[1000];
            sprintf(f1, "visited_temp.txt");
            visited.print_histogram(f1);
            
            sprintf(f1, "ln_dos_temp.txt");
            ln_dos.print_histogram(f1);
            
            if(adjust == true)
            {
                if(translate_particles)
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
                    n_accept_translate = 0;
                    translate_probability  = 0;

                }
                if(cswap_particles)
                {
                     if(cswap_probability < cswap_target_prob && time_current > time_thermo_output)
                     {
                         cswap++;
                         if(cswap > cswap_max)
                             cswap = cswap_max;
                     }
                     else if(cswap_probability > cswap_target_prob && time_current > time_thermo_output)
                     {
                         cswap--;
                         if(cswap < 1)
                         {
                             cswap = 1;
                         }
                     }
                    n_accept_cswap = 0;
                    cswap_probability = 0;

                }
                
                
                
                n_accept_translate = 0;
                n_accept_cswap = 0;
                n_accept_swap = 0;
                n_accept_vol = 0;
                since_last=0;
            }
            
        }
        
        if(time_current%time_xyz_output == 0)
        {
            print_xyz(xyzOut, x, types, L);
        }
        
        if(time_current%time_check_flat == 0 && time_current > time_check_flat)
        {
            double avg_flat=0, min_flat=0, max_flat=0, metric=0;
            bool check_flat =  check_flatness(visited, flatness, &avg_flat, &min_flat, &max_flat, &metric);
            
            std::cout << std::setprecision(5) << "flatness min\t" << min_flat << "\tmax\t" << max_flat << "\tm\t" << metric << std::endl;
            if(check_flat)
            {
                std::cout << "histogram converged with " << flatness << " for f = " << std::setprecision(32) << f << std::endl;
                
            }
            if(check_flat == true)
            {
                is_flat = true;
            }
            
        }
        time_since++;
        time_current++;
        time_local++;
        since_last++;
        
    }
    

    char f1[1000];
    sprintf(f1, "visited_%.10f.txt", f);
    visited.print_histogram(f1);
    sprintf(f1, "ln_dos_%.10f.txt", f);
    ln_dos.print_histogram(f1);

    
    
}

bool WL_2D::run_step(double f_t, bool adjust, bool append_file, bool check_now)
{
    cswap = x.size();

    time_current = 0;
    assert(setup_initialized == true);
    
    std::ofstream xyzOut;
    std::ofstream thermoOut;
    std::ofstream loggerOut;
    
    
    if(append_file)
    {
        xyzOut.open(xyzFilename,std::ios::app);
        thermoOut.open(thermoFilename,std::ios::app);
        loggerOut.open(loggerFilename,std::ios::app);
    }
    else
    {
        xyzOut.open(xyzFilename);
        thermoOut.open(thermoFilename);
        loggerOut.open(loggerFilename);
        
    }
    
    
    f = f_t;
    ln_f = log(f);
    
    ln_dos_old = ln_dos;
    //first calculate the system potential energy
    potential =0;
    calc_pe_brute(x, types, L, &potential, cutoff);
    

    
    int time_since=0;

    bool is_flat = false;
    
    //we need to set up M_current
    if(cswap_particles == true)
        count_type0();
    
    if(vol_change == true)
        calc_vol();
    
    int time_local = 0;
    int n_accept_local;
        
        
    if(translate_particles)
    {
        n_accept_local = translate_WL();
        n_accept_translate += n_accept_local;
    }
    
    nsq_neighbor_check_fast(x, nbr, cutoff, skin, L);
    
    //check to see if we need to rebuild the neighborlist
    if(swap_particles)
    {
        n_accept_local = swap_WL();
        n_accept_swap += n_accept_local;
    }
    
    if(cswap_particles)
    {
        n_accept_local = cswap_WL();
        n_accept_cswap += n_accept_local;
    }
    
    if(vol_change)
    {
        n_accept_local = vol_WL();
        n_accept_vol += n_accept_local;
        
        nsq_neighbor_rebuild(x, nbr, cutoff, skin, L);
    }
    if(check_now == true)
    {
        std::cout << "M_current: " << M_current << std::endl;
        
        
        thermoOut << time_current << "\t" << potential << "\t";
        std::cout << time_current << "\t" << potential << "\t";

        
        double translate_probability = (double)n_accept_translate/(double)(iterations*x.size());
        double cswap_probability = (double)n_accept_cswap/(iterations*cswap);
        double swap_probability = (double)n_accept_swap/(iterations*swap);
        double vol_probability = (double)n_accept_vol/(iterations*vol);
        
        n_accept_translate = 0;
        n_accept_cswap = 0;
        n_accept_swap = 0;
        n_accept_vol =0;
        iterations = 1;
        thermoOut << time_current << "\t" << potential << "\t";
        std::cout << time_current << "\t" << potential << "\t";
        if(translate_particles)
        {
            thermoOut << translate_probability << "\t";
            std::cout << "trans:\t"<< translate_probability << "\t";
        }
        if(swap_particles)
        {
            thermoOut << swap_probability << "\t";
            std::cout << "swap:\t" << swap_probability << "\t";
        }
        if(cswap_particles)
        {
            thermoOut << cswap_probability << "\t";
            std::cout << "cswap:\t" << cswap_probability << "\t";
        }
        if(vol_change)
        {
            thermoOut << vol_probability << "\t";
            std::cout << "vol:\t" << vol_probability << "\t";
        }
        thermoOut << std::endl;
        std::cout << std::endl;
        
        loggerOut << time_current << "\t" << potential << std::endl;
                
        char f1[1000];
        sprintf(f1, "visited_temp.txt");
        visited.print_histogram(f1);
        
        sprintf(f1, "ln_dos_temp.txt");
        ln_dos.print_histogram(f1);
        
        print_xyz(xyzOut, x, types, L);
        double avg_flat, min_flat, max_flat, metric;
        bool check_flat =  check_flatness(visited, flatness, &avg_flat, &min_flat, &max_flat, &metric);
                
        std::cout << std::setprecision(5) << "flatness min\t" << min_flat << "\tmax\t" << max_flat << "\tm\t" << metric << std::endl;
        if(check_flat)
        {
            std::cout << "histogram converged with " << flatness << " for f = " << std::setprecision(32) << f << std::endl;
            
        }
        if(check_flat == true)
        {
            is_flat = true;
            
            char f1[1000];
            sprintf(f1, "visited_%.10f.txt", f);
            visited.print_histogram(f1);
            sprintf(f1, "ln_dos_%.10f.txt", f);
            ln_dos.print_histogram(f1);
        }
        
    }
    iterations++;
    
    return is_flat;
}


//    int translate_WL(coordlist_t& x, std::vector<int>& types, std::vector<neighbor> &nbr, double *L, double dx, double *potential_energy_total, histogram_1D& ln_dos, histogram_1D& visited_local, histogram_1D& ln_dos_local, double ln_factor, double min_energy, double max_energy)

int WL_2D::translate_WL()
{
    
    //precalculate 1/L
    //this could obviously be a global variable, but it's not
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
            double delta_PE = potential_energy_new - potential_energy_old;
            
            double pe_old = potential;
            double pe_new = potential + delta_PE;
            
            //std::cout << delta_PE << "\t" << pe_old << "\t" << pe_new << std::endl;
            
            U_max = ln_dos.get_max(M_current);
            U_min = ln_dos.get_min(M_current);
            

            if(pe_new < U_max  && pe_new > U_min)
            {
                
                double ln_dos_new = ln_dos.get_hist(pe_new, M_current);
                double ln_dos_old = ln_dos.get_hist(pe_old, M_current);
                
                double delta_ln_dos = ln_dos_old-ln_dos_new;
                
                //if the potential energy increases, check to see if we should accept
                double rand_value = drand48();
                if(exp(delta_ln_dos) > rand_value)
                {
                    //I'm just giving it a small weighting factor so it doesn't become really large
                    
                    potential = pe_new;
                    n_accept++;
                    
                    //assign the translated position to the main position array
                    //making sure to apply PBC
                    for( int k=0; k<3; k++)
                    {
                        //new_x[k] -= L[k]*anint(new_x[k]*L_inv[k]);
                        x[i][k] = new_x[k];
                        nbr[i].dx[k] += new_dx[k];
                        
                        
                    }
                    
                    visited.insert(pe_new, M_current);
                    ln_dos.insert(pe_new, M_current, ln_f);
                    
                    
                }
                else
                {
                    visited.insert(pe_old, M_current);
                    ln_dos.insert(pe_old, M_current, ln_f);
                    
                }
            }
            
        }
        
    }
    return n_accept;
    
    
}
int WL_2D::swap_WL()
{
    
    //precalculate 1/L
    double L_inv[3];
    for(int k=0; k<3; k++)
        L_inv[k] = 1.0/L[k];
    
    
    double pot_temp = 0;
    
    int n_accept=0;
    std::vector<int> swapped;
    
    for(int i=0; i<x.size(); i++)
        swapped.push_back(0);
    
    int n_attempts1 = 0;
    int n_attempts2 = 0;
    int max_attempts = 5000;
    int swap_local = swap;
    if(swap > M_current)
    {
        swap_local = M_current;
    }
    for(int i=0; i<swap_local; i++)
    {
        
        bool picking = true;
        
        int n1,t1;
//        while(picking)
//        {
//            
//            n1 = drand48()*x.size();
//            t1 = types[n1];
//            
//            if(swapped[n1] < 1)
//            {
//                picking = false;
//            }
//            n_attempts1++;
//            //we need to give it some way to break from the loop if we get stuck
//            if(n_attempts1 > max_attempts)
//            {
//                //std::cout << n_attempts1 << std::endl;
//                i = swaps+1;
//                break;
//            }
//        }
        
        n1 = drand48()*x.size();
        t1 = types[n1];
            
           
        picking = true;
        int n2, t2;
        while(picking)
        {
            n2 = drand48()*x.size();
            t2 = types[n2];
            
            //std::cout << "n1:\t" << n1 << "\t" << types[n1] << "\t" << swapped[n2] << "\tn2\t: " << n2 << "\t" <<  types[n2] << "\t" << swapped[n2] << "\t i:\t " << i << "\tswap:\t" << swap_local << "\t" << M_current <<  std::endl;
            
            
            if(t1 != t2)
                picking = false;
        
            n_attempts2++;
            //we need to give it some way to break from the loop if we get stuck
            if(n_attempts2 > max_attempts)
            {
                //std::cout << n_attempts2 << std::endl;
                i = swap_local+1;
                break;
            }
        }

//        while(picking)
//        {
//            n2 = drand48()*x.size();
//            t2 = types[n2];
//            
//            std::cout << "n1:\t" << n1 << "\t" << types[n1] << "\t" << swapped[n2] << "\tn2\t: " << n2 << "\t" <<  types[n2] << "\t" << swapped[n2] << "\t i:\t " << i << "\tswap:\t" << swap_local << "\t" << M_current <<  std::endl;
//
//            if(swapped[n2] < 1)
//            {
//                
//                if(t1 != t2)
//                    picking = false;
//            }
//            n_attempts2++;
//            //we need to give it some way to break from the loop if we get stuck
//            if(n_attempts2 > max_attempts)
//            {
//                //std::cout << n_attempts2 << std::endl;
//                i = swaps+1;
//                break;
//            }
//        }
        
        double potential_energy_new = 0.0;
        double potential_energy_old = 0.0;
        bool overlap1 = false;
        bool overlap2 = false;
        
        
        //calculate the potential energy of particle "i" with all of its neighbors "j"
        for(int j=1; j<nbr[n1].member.size(); j++)
        {
            int nbr_id = nbr[n1].member[j];
            
            //calculate the separation between two particles,
            //taking into account periodic boundary conditions
            double rij[3];
            double r2=0;
            double r2_new=0;
            
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
            
            //if one of the particles in the list is the one we are swapping with,
            //we need to swap the type.
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
            double delta_PE = potential_energy_new - potential_energy_old;
            
            double pe_old = potential;
            double pe_new = potential + delta_PE;
            
            //std::cout << delta_PE << "\t" << pe_old << "\t" << pe_new << std::endl;
            U_max = ln_dos.get_max(M_current);
            U_min = ln_dos.get_min(M_current);
            
            if(pe_new < U_max  && pe_new > U_min)
            {
                
                double ln_dos_new = ln_dos.get_hist(pe_new, M_current);
                double ln_dos_old = ln_dos.get_hist(pe_old, M_current);
                
                double delta_ln_dos = ln_dos_old-ln_dos_new;
                
                
                
                //if the potential energy increases, check to see if we should accept
                double rand_value = drand48();
                if(exp(delta_ln_dos) > rand_value)
                {
                    //I'm just giving it a small weighting factor so it doesn't become really large
                    
                    potential = pe_new;
                    n_accept++;
                    
                    types[n1] = t2;
                    types[n2] = t1;
                    
                    
                    visited.insert(pe_new, M_current);
                    ln_dos.insert(pe_new, M_current, ln_f);
                    
                    swapped[n1] = 1;
                    swapped[n2] = 1;
                }
                else
                {
                    
                    visited.insert(pe_old, M_current);
                    ln_dos.insert(pe_old, M_current, ln_f);
                    
                }
            }
        }
        
    }
    return n_accept;
    
    
}
int WL_2D::count_type0()
{
    //in this case M_current is the fraction of type 0
    int type0=0;
    for(int i=0; i<types.size(); i++)
    {
        if(types[i] == 0)
            type0++;
    }
    M_current = (double)type0;
    return type0;
}


//composition swap move
int WL_2D::cswap_WL()
{
    //in this case M_current is the fraction of type 0
    int type0 = count_type0();
    
    //precalculate 1/L
    double L_inv[3];
    for(int k=0; k<3; k++)
        L_inv[k] = 1.0/L[k];
    
    
    double pot_temp = 0;
    
    int n_accept=0;
    std::vector<int> swapped;
    
    for(int i=0; i<x.size(); i++)
        swapped.push_back(0);
    int i=0;
    while(i < cswap)
    {
        int type0_temp = M_current;
        int n1,t1;
        
        n1 = x.size()*drand48();
        if(swapped[n1] == 0)
        {
            i++;
            t1 = types[n1];
      
           
            int t2 = swap_type(t1);
            
            
            if(t1 == 0)
                type0_temp--;
            else
                type0_temp++;
            
            
            
            double M_current_new = (double)type0_temp;
            
            double M_current_old = M_current;
            double potential_energy_new = 0.0;
            double potential_energy_old = 0.0;
            bool overlap1 = false;
            
            
            //calculate the potential energy of particle "i" with all of its neighbors "j"
            for(int j=1; j<nbr[n1].member.size(); j++)
            {
                int nbr_id = nbr[n1].member[j];
                
                //calculate the separation between two particles,
                //taking into account periodic boundary conditions
                double rij[3];
                double r2=0;
                double r2_new=0;
                
                for(int k=0; k<3; k++)
                {
                    rij[k] = x[n1][k] - x[nbr_id][k];
                    
                    double pbc  = L[k]*anint(rij[k]*L_inv[k]);
                    
                    rij[k] -= pbc;
                    r2 += rij[k]*rij[k];
                }
                
                pot_temp = 0;
                //assume no overlaps in this configuration
                //we'll need to modify this to do two different sized particles, if we want to do that
                LJ_potential(r2, t1, types[nbr_id], &pot_temp);
                potential_energy_old += pot_temp;
                
                int tt = types[nbr_id];
                
                
                overlap1 = LJ_potential(r2, t2, tt, &pot_temp);
                if(overlap1==true)
                {
                    break;
                }
                potential_energy_new += pot_temp;
                
            }
            if(overlap1 == false)
            {
                double delta_PE = potential_energy_new - potential_energy_old;
                
                double pe_old = potential;
                double pe_new = potential + delta_PE;
                
                //std::cout << delta_PE << "\t" << pe_old << "\t" << pe_new << std::endl;
                if(M_current_new < M_max && M_current_new >= M_min)
                {
                    U_max = ln_dos.get_max(M_current_new);
                    U_min = ln_dos.get_min(M_current_new);
                
                    if(pe_new < U_max  && pe_new > U_min)
                    {
                        
                        double ln_dos_new = ln_dos.get_hist(pe_new, M_current_new);
                        double ln_dos_old = ln_dos.get_hist(pe_old, M_current_old);
                        
                        double delta_ln_dos = ln_dos_old-ln_dos_new;
                        
                        
                        
                        //if the potential energy increases, check to see if we should accept
                        double rand_value = drand48();
                        if(exp(delta_ln_dos) > rand_value)
                        {
                            //I'm just giving it a small weighting factor so it doesn't become really large
                            
                            potential = pe_new;
                            n_accept++;
                            
                            types[n1] = t2;
                            
                            
                            visited.insert(pe_new, M_current_new);
                            ln_dos.insert(pe_new, M_current_new, ln_f);
                            
                            swapped[n1] = 1;
                            
                            type0 = count_type0();

                        }
                        else
                        {
                            
                            visited.insert(pe_old, M_current_old);
                            ln_dos.insert(pe_old, M_current_old, ln_f);
                            
                        }
                    }
                }
            }
        }
    }
    return n_accept;
    
    
}
double WL_2D::calc_vol()
{
    M_current = log(L[0]*L[1]*L[2]);
    return M_current;
}



int WL_2D::vol_WL()
{
 
    //this will calculate the volume and set M_current to this volume
    
    int n_accept  = 0;
    coord_t L_new;
    for(int k=0; k<3; k++)
        L_new.push_back(L[k]);

    //vol corresponds to the number of attempted volume changes per cycle
    //in WL we could probably do a bunch of them each cycle
    //since ultimately we just want to sample different configurations
    for(int i=0; i<vol; i++)
    {
        double vol_old = L[0]*L[1]*L[2];
        double M_current_old = calc_vol();
        
        double pe_old = potential;
        double pe_new;
        
        double lnvolnew = log(vol_old)+(2.0*drand48()-1.0)*ln_vol_change_max;
        double vol_new = exp(lnvolnew);
        
        double scaling_factor = pow(vol_new/vol_old, 1.0/3.0);
        
        coordlist_t x_new = x;
        
        for(int j=0; j<x.size(); j++)
        {
            for(int k=0; k<3; k++)
                x_new[j][k] *= scaling_factor;
        }
        for(int k=0; k<3; k++)
            L_new[k] = scaling_factor*L[k];
        
        bool overlap = calc_pe(x_new, types, L_new, &pe_new, cutoff);
        

        if(overlap == false)
        {
            double M_current_new = lnvolnew; //vol_new;
            
            if(M_current_new < M_max && M_current_new >= M_min)
            {
                U_max = ln_dos.get_max(M_current_new);
                U_min = ln_dos.get_min(M_current_new);
                
                if(pe_new < U_max  && pe_new > U_min)
                {
                    
                    double ln_dos_new = ln_dos.get_hist(pe_new, M_current_new);
                    double ln_dos_old = ln_dos.get_hist(pe_old, M_current_old);
                    
                    double delta_ln_dos = (ln_dos_old-log(vol_old)) - (ln_dos_new-log(vol_new));
                    
                    
                    
                    //if the potential energy increases, check to see if we should accept
                    double rand_value = drand48();
                    if(exp(delta_ln_dos) > rand_value)
                    {
                        //I'm just giving it a small weighting factor so it doesn't become really large
                        
                        potential = pe_new;
                        n_accept++;
                        
                        x = x_new;
                        
                        L[0] = L_new[0];
                        L[1] = L_new[1];
                        L[2] = L_new[2];
                        
                        visited.insert(pe_new, M_current_new);
                        ln_dos.insert(pe_new, M_current_new, ln_f);
                        
                        
                    }
                    else
                    {
                        
                        visited.insert(pe_old, M_current_old);
                        ln_dos.insert(pe_old, M_current_old, ln_f);
                        
                    }
                }
            }

        }
        
        
    }
    return n_accept;
}
