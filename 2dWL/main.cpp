/*
 *  main.cpp
 *  nve_LJ
 *
 *  Created by Christopher Iacovella on 8/1/11.
 *  A bare bones code that performs an MC simulation of LJ particles

 */

#include "main.h"

//this routine will first run a bunch of 1d WL to figure out the minimum accessible energy

void calc_min_energy(int argc , char **argv)
{
    
    MPI_Init(&argc,&argv);
    //MPI::Init(argc,argv);

    double cutoff = 2.5;
    double skin = 1.0;
    
    //MC parameters
	double dx = skin/3.1;
    double dx_max = skin/2.0001; // max displacement
    double dx_min = 0.001;
	double dx_target = 0.5;	//target acceptance probability

    
	//System parameters
	int N_particles = 128;
    int Na_particles_min = 120; 
    int Na_particles = N_particles;
    int N_bin = 1;
        
    
	double number_density = 0.5;
    
    
    //we need to define the maximum
    double U_max = -200.0;
    double U_min = -700.0;
    double U_min2 = -600.0;
    double U_bin = 2.0;

    int time_check_flat = 1000;
    int max_loop = (N_particles-Na_particles_min)/N_bin;
    
    
    std::ofstream minOut("min_energy.txt");
    
    
    for(int i=0; i<max_loop; i+=1)
    {
        //MPI::COMM_WORLD.Barrier();
        MPI_Barrier(MPI_COMM_WORLD);
        
        //we'll first run each system a short while in standard NVT to get a configuration
        //in the desired energy range
        metropolis_NVT mNVT;
        
        
        int threads = mNVT.get_threads();
        int rank = mNVT.get_rank();
        
        short unsigned int sseed = 12345+rank*100+rank*11+rank;
        //srand48(sseed);
        //let's just have every system start from the same configuration and then diverge
        
        int Na_current = N_particles-i*N_bin;
        std::cout << rank << "\t" << Na_current << std::endl;

        
        //define arrays to hold position (x)
        //note coordlist_t is just a convenient container class I defined that is essentially a double array x[n_particles][3]
        coordlist_t x;
        std::vector<int> types;
        
        //define our box array
        double L[3];
        init_system_binary(x, N_particles, Na_current, number_density, L, types);
        std::cout << "system with " << x.size() << " particles initialized, Na: "<< Na_current << std::endl;
    
        
        double T_min = 0.1;
        double T_max = 6.0;
        
        double T = T_min+((double)rank/(double)threads)*(T_max-T_min); //just a function to give us different starting temperatures for each walker
        
        mNVT.init_translate(dx, dx_min, dx_max, dx_target);
        mNVT.init_system(L,  sseed, T);
        mNVT.set_configuration(x, types);
        mNVT.init_nlist(cutoff, skin);
        mNVT.setup();
        
        
        
        int time_U_update = 100;
        int time_xyz_output = 10000;
        int time_thermo_output = 5000;
        
        
        
        char f1[1000];
        sprintf(f1, "traj_NVT_%d.xyz", rank);
        std::ofstream xyzOut(f1);
        
        sprintf(f1, "thermo_NVT_%d.txt", rank);
        std::ofstream thermoOut(f1);
        
        int time_mc = 5000;		//total number of mc timesteps for the
        
        mNVT.converge_window(time_mc, time_U_update, time_xyz_output, xyzOut, time_thermo_output, thermoOut, true, U_min2, U_max, T_min, T_max);
        mNVT.get_configuration(x, types);
        
        //MPI::COMM_WORLD.Barrier();
        MPI_Barrier(MPI_COMM_WORLD);

        if(rank == 0)
        {
            std::cout << "**************" << std::endl;
            std::cout << "initial setup done, starting WL" << std::endl;
            std::cout << "**************" << std::endl;
        }
        std::cout << "rank: " << rank << "\tpot:\t" << mNVT.potential << std::endl;
    
        //MPI::COMM_WORLD.Barrier();
        MPI_Barrier(MPI_COMM_WORLD);

   
        
        double f = exp(1);

        double flatness = 0.5;
        WL_NVT_1D WL;
        WL.init_translate(dx, dx_min, dx_max, dx_target);
        //if(i > 0)
        //    WL.init_swap(10, N_particles/2.0, 0.05);
        WL.init_system(L,  sseed);
        WL.init_nlist(cutoff, skin);
        WL.init_WL(U_min, U_max, U_bin, flatness);
        WL.set_configuration(x, types);
        WL.setup();
        
        
        sprintf(f1, "traj_WL_%d.xyz", rank);
        std::ofstream xyzOut2(f1);
        
        sprintf(f1, "thermo_WL_%d.txt", rank);
        std::ofstream thermoOut2(f1);
        
        WL.flatness = flatness;
        WL.set_max_time(500000);
        

        WL.run(f, time_U_update, time_check_flat, time_xyz_output, xyzOut2, time_thermo_output, thermoOut2, true);
        

        if(rank ==0)
        {
            std::cout << "=============" << std::endl;
            std::cout << "Na:\t" << Na_current << "\tpe_min:\t" << WL.min_energy() << std::endl;
            minOut << Na_current << "\t" << WL.min_energy() << std::endl;
            std::cout << "=============" << std::endl;
            
        }
        //MPI::COMM_WORLD.Barrier();
        MPI_Barrier(MPI_COMM_WORLD);

    }
    
   	MPI_Finalize ( );
    
}

void WL_2D(int argc , char **argv)
{
    //needed for MPI
    MPI_Init(&argc,&argv);
    
    //the first step in this code will be to generate configurations within the desired energy range
    //to do this, we'll first run each system in standard NVT metropolis MC
    //adjusting the temperature to hit the target energy
    
    metropolis_NVT mNVT;
    
    //calculate the total number of threads and the current process rank
    //this can obviously be called directly via MPI
    int threads = mNVT.get_threads();
    int rank = mNVT.get_rank();
    
    //we'll use processor rank to get each thread a unique seed
    short unsigned int sseed = 12345+rank*100+rank*11+rank;
    srand48(sseed);
    
    
    //===================
	//System parameters
    //===================
    
    
    //pair potential cutoff
    double cutoff = 2.5;
    
    //the code uses a brute force neighborlist for efficiency
    //the skin can be set rather large, even if cutoff+skin > half the box length
    //since a particle can only be added to the neighborlist once
    //however, the max translation is limited by the skin size
    double skin = 1.0;
    
    //MC parameters
	double dx = skin/3.1;           // our initial dx, we will adjust it later
    double dx_max = skin/2.0001;    // max displacement
    double dx_min = 0.001;          // min displacement
	double dx_target = 0.5;         //target acceptance probability
    
    
    //system size
    int N_particles = 128;
    
    //we also need to define the minimum number of type A particles
    int Na_particles_min = 64;
    int Na_particles_max = 128;

    //we have a binary mixtures, so let's define how many are type A
    int Na_particles = Na_particles_min+(Na_particles_max-Na_particles_min)*rank/threads;

    //a quick output to make sure things are as we expect
    std::cout << rank << "\t" << Na_particles << std::endl;
    
    //number density of the system
	double number_density = 0.75;
    
    //define arrays to hold position and type
	coordlist_t x;
    std::vector<int> types;
    
	//define our box array
	double L[3];
    
    //randomly insert particles in a box
    //if dense systems are needed, we may have to add a routine to compress the box
    init_system_binary(x, N_particles, Na_particles, number_density, L, types);
    
    //===================
	//MC parameters
    //===================
    
    //we need to define an initial range to converge the system
    //this may or may not be the overall energy range for the WL routine,
    //i.e., U_min may need to be a function of the 2nd dimension
    double U_min = -200.0;
    double U_max = -100.0;
    
    //this is the temperature range we will explore in order to make sure we are in the energy range
    double T_min = 0.1;
    double T_max = 6.0;
    
    //just a function to give us different starting temperatures for each walker
    double T = T_min+((double)rank/(double)threads)*(T_max-T_min); 
    
    //for efficiency, we calcuate system energy by looking at energy differences
    //to avoid round-off issues, we periodically recalculate the energy via brute force calculation
    int time_U_update = 100;
    
    //these are just time parameters for outputting data
    int time_xyz_output = 10000;
    int time_thermo_output = 5000;
    
    //output files for MC
    char f1[1000];
    sprintf(f1, "traj_NVT_%d.xyz", rank);
    std::ofstream xyzOut(f1);
    
    sprintf(f1, "thermo_NVT_%d.txt", rank);
    std::ofstream thermoOut(f1);
    
    int time_mc = 5000;		//total number of mc timesteps for the
    

    
    //setup MC parameters
    //note, if you don't define everything that is needed
    //the code should alert you and terminate
    mNVT.init_translate(dx, dx_min, dx_max, dx_target);
    mNVT.init_system(L,  sseed, T);
    mNVT.set_configuration(x, types);
    mNVT.init_nlist(cutoff, skin);
    mNVT.setup();
    
    //run the MC routine to converge to the energy window
    mNVT.converge_window(time_mc, time_U_update, time_xyz_output, xyzOut, time_thermo_output, thermoOut, true, U_min, U_max, T_min, T_max);
    
    //update the configuration and types
    mNVT.get_configuration(x, types);

    
    MPI_Barrier(MPI_COMM_WORLD);

    
    if(rank == 0)
    {
        std::cout << "**************" << std::endl;
        std::cout << "initial setup done, starting WL" << std::endl;
        std::cout << "**************" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    
    //===================
	//WL parameters
    //===================
    
    
    //flatness criteria
    double flatness = 0.75;
    
    //how frequently should we check to see if the histogram is flat
    int time_check_flat = 1000;
    
    //how often to output thermo
    time_thermo_output = 1000;
    
    //thermo and configuration coordinate file streams
    sprintf(f1, "traj_WL_%d.xyz", rank);
    std::ofstream xyzOut2(f1);
    
    sprintf(f1, "thermo_WL_%d.txt", rank);
    std::ofstream thermoOut2(f1);
    
    //initial modification factor
    double f = exp(1);
    
    //set up energy range and bin size
    U_min = -600.0;
    U_max = -100.0;
    double U_bin = 2.0;
    
    //since this is 2D, we need to define the range
 
    double N_min = Na_particles_min;
    double N_max = Na_particles_max+1;   //N_max always needs +1 to ensure we include the upper window
    double N_bin = 1;
 
    
    //===================
	//WL Setup
    //===================
    
    WL_NVT_2D WL;
    WL.init_translate(dx, dx_min, dx_max, dx_target);
    WL.init_cswap(1, N_particles, 0.05);
    WL.init_system(L,  sseed);
    WL.init_nlist(cutoff, skin);
    
    //we first initialize uniform, then will run for a while to calculate U min
    WL.init_uniform_WL(U_min, U_max, U_bin, flatness, N_min, N_max, N_bin);
    WL.set_configuration(x, types);
    WL.setup();
    
    int WL_run_time = 100000;
    WL.calc_U_min(f, WL_run_time, time_U_update, time_check_flat, time_thermo_output, true);
    
    
    while(f > exp(1e-8))
    {
        
        WL.run(f, time_U_update, time_check_flat, time_xyz_output, xyzOut2, time_thermo_output, thermoOut2, true);
        f = sqrt(f);
        WL.reset_configuration();
        MPI_Barrier(MPI_COMM_WORLD);

    }
    //for the last stage, we won't update dx and we will merge lndos every step to ensure better consistency
    WL.run(f, time_U_update, time_check_flat, time_xyz_output, xyzOut2, time_thermo_output, thermoOut2, false);
    MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize ( );
    
}



void WL_2D_restart(int argc , char **argv)
{
    //needed for MPI
    //MPI::Init(argc,argv);
    MPI_Init(&argc,&argv);

    //we'll restart from the xyz files we output earlier
    //so we don't need to waste time running MC again.
    
    WL_NVT_2D WL;
    
    //calculate the total number of threads and the current process rank
    //this can obviously be called directly via MPI
    int threads = WL.get_threads();
    int rank = WL.get_rank();
    
    
    //read in files from MC
    char f1[1000];
    sprintf(f1, "traj_NVT_%d.xyz", rank);

    
    
    //system size
    int N_particles = 128;
    
    //we also need to define the minimum number of type A particles
    int Na_particles_min = 64;
    int Na_particles_max = 128;
    
    
    //define arrays to hold position and type
	coordlist_t x;
    std::vector<int> types;
    
    load_xyz(f1, x, types);
    
    WL.set_configuration(x, types);

    
    //number density of the system
	double number_density = 0.75;
    
	//define our box array
	double L[3];
    
    double L_temp = pow((N_particles/number_density), 1.0/3.0);
    L[0] = L[1] = L[2] = L_temp;
    
    //we'll use processor rank to get each thread a unique seed
    short unsigned int sseed = 12345+rank*100+rank*11+rank;
    srand48(sseed);
    
    
    
    
    int Na_particles  = WL.count_type0();
    std::cout << rank << "\t" << x.size() << "\t" << Na_particles << std::endl;
    //===================
	//System parameters
    //===================
    
    //pair potential cutoff
    double cutoff = 2.5;
    
    //the code uses a brute force neighborlist for efficiency
    //the skin can be set rather large, even if cutoff+skin > half the box length
    //since a particle can only be added to the neighborlist once
    //however, the max translation is limited by the skin size
    double skin = 1.0;
    
    //MC parameters
	double dx = skin/3.1;           // our initial dx, we will adjust it later
    double dx_max = skin/2.0001;    // max displacement
    double dx_min = 0.001;          // min displacement
	double dx_target = 0.5;         //target acceptance probability
    
    
    //for efficiency, we calcuate system energy by looking at energy differences
    //to avoid round-off issues, we periodically recalculate the energy via brute force calculation
    int time_U_update = 100;
    
    //these are just time parameters for outputting data
    int time_xyz_output = 10000;
    int time_thermo_output = 5000;
    

    
    int time_mc = 5000;		//total number of mc timesteps for the
    

    WL.init_ln_dos_from_file("ln_dos_1.0000076294.txt");
    
    MPI_Barrier(MPI_COMM_WORLD);

    
    //===================
	//WL parameters
    //===================
    
   /*
    //flatness criteria
    double flatness = 0.75;
    
    //how frequently should we check to see if the histogram is flat
    int time_check_flat = 1000;
    
    //how often to output thermo
    time_thermo_output = 1000;
    
    //thermo and configuration coordinate file streams
    sprintf(f1, "traj_WL_%d.xyz", rank);
    std::ofstream xyzOut2(f1);
    
    sprintf(f1, "thermo_WL_%d.txt", rank);
    std::ofstream thermoOut2(f1);
    
    //initial modification factor
    double f = exp(1);
    
    //set up energy range and bin size
    U_min = -600.0;
    U_max = -100.0;
    double U_bin = 2.0;
    
    //since this is 2D, we need to define the range
    
    double N_min = Na_particles_min;
    double N_max = Na_particles_max+1;   //N_max always needs +1 to ensure we include the upper window
    double N_bin = 1;
    
    
    //===================
	//WL Setup
    //===================
    
    WL.init_translate(dx, dx_min, dx_max, dx_target);
    WL.init_cswap(1, N_particles, 0.05);
    WL.init_system(L,  sseed);
    WL.init_nlist(cutoff, skin);
    
    //we first initialize uniform, then will run for a while to calculate U min
    WL.init_uniform_WL(U_min, U_max, U_bin, flatness, N_min, N_max, N_bin);
    WL.set_configuration(x, types);
    WL.setup();
    
    int WL_run_time = 100000;
    WL.calc_U_min(f, WL_run_time, time_U_update, time_check_flat, time_thermo_output, true);
    
    
    while(f > exp(1e-8))
    {
        
        WL.run(f, time_U_update, time_check_flat, time_xyz_output, xyzOut2, time_thermo_output, thermoOut2, true);
        f = sqrt(f);
        WL.reset_configuration();
        MPI::COMM_WORLD.Barrier();
        
    }
    //for the last stage, we won't update dx and we will merge lndos every step to ensure better consistency
    WL.run(f, time_U_update, time_check_flat, time_xyz_output, xyzOut2, time_thermo_output, thermoOut2, false);
    MPI::COMM_WORLD.Barrier();
    */
	MPI_Finalize ( );
    
}



int main(int argc , char **argv)
{

	//calc_min_energy(argc , argv);
    //WL_2D_restart(argc , argv);

    WL_2D(argc , argv);
	return 0;
}
