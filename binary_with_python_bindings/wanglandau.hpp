/*
 *  wanglandau.hpp
 *  simple_LJ
 *
 *  Created by Christopher Iacovella
 *  Copyright 2012. All rights reserved.
 *
 */

#ifndef _wanglandau_hpp
#define _wanglandau_hpp

class WL_2D{
    
    
public:
    //translation parameters
    double dx;
    double dx_max;
    double dx_min;
    double dx_target_prob;
    bool translate_particles;
    int nn;
    
    //identiy swap parameters
    int swap;
    int swap_max;
    double swap_target_prob;
    bool swap_particles;
    
    
    //composition swap parameters
    int cswap; //initial number of times we will try to swap
    int cswap_max; //maximum times we will try to swap
    double cswap_target_prob; //target probability for swaps
    bool cswap_particles;
    
    //volume swap parameters
    int vol;    //number of volume change attempts per cycle
    int vol_max;    //max number of volume change attempts per cycle
    
    double ln_vol_change_max;
    double vol_target_prob;
    bool vol_change;
    
    
    //system parameters
    coord_t L;
    short unsigned int sseed;
    int time_current;
    int rank;
    int threads;
    double potential;
    bool system_initialized;
    
    //configuration
    int N_particles;
    coordlist_t x;
    std::vector<int> types;
    
    coordlist_t x_init;
    types_t types_init;
    coord_t L_init;
    
    bool configuration_initialized;
    
    //nlist/interaction parameters
    double skin;
    double cutoff;
    
    std::vector<neighbor> nbr;
    
    bool nlist_initialized;
    
    //we'll have generic variables that apply to the 2nd dimension
    //so we can use the same code for volume changes or composition changes
    //let's call this variable M
    double M_current;
    double M_min;
    double M_max;
    double M_bin;
    
    
    //1D WL parameters
    double U_min;
    double U_max;
    double U_max_global;
    double U_bin;
    double U_bin_global;
    

    
    double flatness;
    double f;
    double ln_f;
    histogram_2D visited;
    histogram_2D ln_dos;
    histogram_2D ln_dos_old;
    
    bool WL_initialized;
    
    
    //file output
    char xyzFilename[1000];
    char thermoFilename[1000];
    char loggerFilename[1000];

    
    bool xyzOut_initialized;
    bool thermoOut_initialized;
    bool loggerOut_initialized;

    void set_thermo_filename(const char* filename);
    void set_traj_filename(const char* filename);
    
    bool setup_initialized;
    
    
    int n_accept_translate;
    int n_accept_swap;
    int n_accept_cswap;
    int n_accept_vol;
    int iterations;
    int since_last;
    WL_2D();
    
    void init_uniform_WL(double U_min_t, double U_max_t, double U_bin_t, double flatness_t, double M_min_t, double M_max_t, double M_bin_t);
    
    void init_from_params_WL(const char* filename, double flatness_t);
    
    void init_ln_dos_from_file(const char* filename_ln_dos);
    
    void init_y(double M_min_t, double M_max_t, double M_bin_t);
    
    void push_back_x(double U_min_t, double U_max_t, double U_bin_t);
    void setup_histograms(double flatness_t);
    
    void init_translate(double dx_t, double dx_min_t, double dx_max_t, double dx_target_prob_t);
    
    void init_system(coord_t L_t,  short unsigned int sseed_t);
    
    void init_swap(int swap_t, int swap_max_t, double swap_target_prob_t);
    void init_cswap(int cswap_t, int cswap_max_t, double cswap_target_prob_t);;
    
    void init_vol_change(int vol_t, int vol_max_t, double vol_target_prob_t, double ln_vol_change_max_t);
    
    
    void init_nlist(double cutoff_t, double skin_t);
    void set_configuration(coordlist_t x_t, types_t types_array_t);
    void reset_configuration();
    
    
    
    void get_configuration(coordlist_t& x_t, types_t& types_array_t);
    void setup();

    void calc_U_min(double f_t, int run_time, int time_U_update, int time_check_flat, int time_thermo_output, bool adjust, bool append_file);
    bool run_step(double f_t, bool adjust, bool append_file, bool check_now);
    void run(double f_t, int time_U_update, int time_check_flat, int time_xyz_output, int time_thermo_output, bool adjust, bool append_file);
    void clear_visited();
    
    //translate move
    int translate_WL();
    
    //identifty swap
    int swap_WL();
    
    //composition swap move
    int count_type0();
    int cswap_WL();
    
    //volume change
    double calc_vol();
    int vol_WL();
  
};
#endif
