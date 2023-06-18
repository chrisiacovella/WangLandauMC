/*
 *  metropolis.hpp
 *  simple_LJ
 *
 *  Created by Christopher Iacovella on 1/25/12.
 *  Copyright 2012. All rights reserved.
 *
 */


//This is a class to handle standard metropolis MC
//This class will contain all the basic routines to minimize the amount of "stuff" in the main function


#ifndef _metropolis_hpp
#define _metropolis_hpp
class metropolis_NVT{
  
    
public:
    //translation parameters
    double dx;
    double dx_max;
    double dx_min;
    double dx_target_prob;
    bool translate_particles;
    
    
    //swap parameters
    int swap;
    double swap_target_prob;
    bool swap_particles;
    
    
    //system parameters
    coord_t L;
    short unsigned int sseed;
    double T;
    double beta;
    int time_current;
    double potential;
    
    bool system_initialized;
    
    //configuration
    int N_particles;

    coordlist_t x;
    types_t types;
    
    bool configuration_initialized;

    //nlist/interaction parameters
    double skin;
    double cutoff;

    std::vector<neighbor> nbr;
    
    bool nlist_initialized;
    bool setup_initialized;

    //file output
    char xyzFilename[1000];
    char thermoFilename[1000];
    
    bool xyzOut_initialized;
    bool thermoOut_initialized;
   
    metropolis_NVT();
    void init_translate(double dx_t, double dx_min_t, double dx_max_t, double dx_target_prob_t);
    void init_system(coord_t L_t,  short unsigned int sseed_t, double T_t);
    void init_swap(int swap_t, double swap_target_prob_t);
    void init_nlist(double cutoff_t, double skin_t);
    void set_configuration(coordlist_t x_t, types_t types_array_t);
    void get_configuration(coordlist_t& x_t, types_t& types_array_t);
    void set_traj_filename(const char* filename);
    void set_thermo_filename(const char* filename);
    void set_temperature(double T_t);   //an additional routine to reset the temperature in case we are cooling or something

    void setup();
    void converge_window(int time_run, int time_U_update, int time_xyz_output, int time_thermo_output, bool adjust, double U_min, double U_max, double T_min, double T_max);
    void run(int time_run, int time_U_update, int time_xyz_output, int time_thermo_output, bool adjust, bool append_file);
    int swap_NVT(coordlist_t& x, std::vector<int>& types, std::vector<neighbor> &nbr, coord_t L, int n_swaps, double *potential_energy_total, double beta);
    int translate_NVT();
};

#endif
