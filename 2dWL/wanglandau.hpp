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

class WL_NVT_1D{
  
    
public:
    //translation parameters
    double dx;
    double dx_max;
    double dx_min;
    double dx_target_prob;
    bool translate_particles;
    int nn;
    
    //swap parameters
    int swap;
    int swap_max;
    double swap_target_prob;
    bool swap_particles;
    
    
    //system parameters
    double L[3];
    short unsigned int sseed;
    int time_current;
    int rank;
    int threads;
    double potential;
    bool system_initialized;
    bool check_max_time;
    int max_time;
    
    //configuration
    int N_particles;
    coordlist_t x;
    std::vector<int> types;
    
    coordlist_t x_init;
    std::vector<int> types_init;
    
    bool configuration_initialized;

    //nlist/interaction parameters
    double skin;
    double cutoff;

    std::vector<neighbor> nbr;
    
    bool nlist_initialized;

    //WL parameters
    double U_min;
    double U_max;
    double U_bin;
    
    
    double flatness;
    double f;
    double ln_f;
    histogram_1D visited;
    histogram_1D ln_dos;
    histogram_1D ln_dos_old;
    
    bool WL_initialized;

    std::ofstream loggerOut;
    
    bool setup_initialized;
    
    double pe_min;
    
    WL_NVT_1D()
    {
        translate_particles = false;
        system_initialized = false;
        nlist_initialized = false;
        configuration_initialized = false;
        swap_particles = false;
        WL_initialized = false;
        setup_initialized = false;
        nn=0;
        time_current = 0;
        pe_min = 1000;
        check_max_time = false;
    }

    void init_WL(double U_min_t, double U_max_t, double U_bin_t, double flatness_t)
    {
        U_min = U_min_t;
        U_max = U_max_t;
        U_bin = U_bin_t;
        flatness = flatness_t;
        
        
        visited.init(U_min, U_max, U_bin);
        
        ln_dos.init(U_min, U_max, U_bin);

        WL_initialized = true;

    }
    void init_translate(double dx_t, double dx_min_t, double dx_max_t, double dx_target_prob_t)
    {
        dx = dx_t;
        dx_max = dx_max_t;
        dx_min = dx_min_t;
        dx_target_prob = dx_target_prob_t;
        
        translate_particles = true;
    }
    
    void init_system(double *L_t,  short unsigned int sseed_t)
    {
        for(int k=0; k<3; k++)
            L[k] = L_t[k];
        
        sseed = sseed_t;
        system_initialized = true;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        MPI_Comm_size(MPI_COMM_WORLD,&threads);
    }
    
    void init_swap(int swap_t, int swap_max_t, double swap_target_prob_t)
    {
        swap = swap_t;
        swap_max = swap_max_t;
        swap_target_prob = swap_target_prob_t;
        swap_particles = true;

        
    }
    void init_nlist(double cutoff_t, double skin_t)
    {
        cutoff = cutoff_t;
        skin = skin_t;
        
        nlist_initialized = true;

    }
    void set_configuration(coordlist_t x_t, std::vector<int> types_t)
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

                types.push_back(types_t[i]);
                types_init.push_back(types_t[i]);
                
            }
        }
        else
        {
            for(int i=0; i< x_t.size(); i++)
            {
                x[i] = x_t[i];
                x_init[i] = x_t[i];
                types[i] = types_t[i];
                types_init[i] = types_t[i];
            }
        }
        configuration_initialized = true;
    }
    void reset_configuration()
    {
        assert(x.size() == x_init.size() );
        
        for(int i=0; i<x.size(); i++)
        {
            x[i] = x_init[i];
            types[i] = types_init[i];
        }
        nsq_neighbor_rebuild(x, nbr, cutoff, skin, L);
        
        
        
    }
    
    
    
    void get_configuration(coordlist_t& x_t, std::vector<int>& types_t)
    {
       if(x.size() != x_t.size())
       {
            x_t.clear();
           types_t.clear();
       
           for(int i=0; i< x.size(); i++)
           {
               x_t.push_back(x[i]);
               types_t.push_back(types[i]);
           }
       }
       else
       {
           for(int i=0; i< x_t.size(); i++)
           {
               x_t[i] = x[i];
               types_t[i] = types[i];
            }
        }
    }
    void setup()
    {
        assert(dx < skin/2.0);
        for(int k=0; k<3; k++)
            assert(cutoff <= L[k]/2.0);
        assert(system_initialized);
        assert(nlist_initialized);
        assert(configuration_initialized);
        assert(WL_initialized);
        assert(x.size() > 0);
        nsq_neighbor_init(x, nbr, cutoff, skin, L);

        potential =0;
        calc_pe_brute(x, types, L, &potential, cutoff);
        if(translate_particles)
            nn+=x.size();
        
        char f1[1000];
        sprintf(f1, "pe_rank%d.txt", rank);
        loggerOut.open(f1);
        
        setup_initialized = true;
        
        
    }
    int get_threads()
    {
        MPI_Comm_size(MPI_COMM_WORLD,&threads);
        return threads;
    }
    int get_rank()
    {
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        return rank;
    }
    
    void init_ln_dos_from_file(const char* filename_ln_dos)
    {
        if(WL_initialized == true)
        {
            ln_dos.read_from_file(filename_ln_dos, get_threads());
        }
        else
        {
            std::cerr<< "please define energy range and binsize first in init_WL" << std::endl;
            assert(WL_initialized==true);
        }
    }

    
    
    void reduce_dos()
    {

        
        visited.splat(0);
        ln_dos.splat(0);
        MPI_Allreduce(visited.array_local,visited.array,visited.size(),MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(ln_dos.array_local,ln_dos.array,ln_dos.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        //MPI::COMM_WORLD.Allreduce(visited.array_local,visited.array,visited.size(),MPI::DOUBLE,MPI::SUM);
        //MPI::COMM_WORLD.Allreduce(ln_dos.array_local,ln_dos.array,ln_dos.size(),MPI::DOUBLE,MPI::SUM);



    }
    void check_reduce()
    {
        double sum_S = 0;
        double sum_V = 0;

        for(int i=0; i<ln_dos.size(); i++)
        {
            sum_S+=ln_dos.array[i]/ln_dos.size();
            sum_V+=visited.array[i]/visited.size();
        }
        std::cout <<"check reduce:\t" << rank << "\t" << sum_S << "\t" << sum_V << std::endl;
    }
    void set_max_time(int max_time_t)
    {
        max_time = max_time_t;
        check_max_time = true;
    }
    
    double min_energy()
    {
        
        for(int i=0; i<visited.size(); i++)
        {
            if(visited.array[i] >0 )
            {
                double val = ln_dos.get_value(i);
                return val;
            }
        }
        return U_max;
    }
       

    
    void run(double f_t, int time_U_update, int time_check_flat, int time_xyz_output, std::ostream& xyzOut, int time_thermo_output, std::ostream& thermoOut, bool adjust)
    {
        assert(setup_initialized == true);
        
        f = f_t;
        ln_f = log(f);

        ln_dos_old = ln_dos;
        //first calculate the system potential energy
        potential =0;
        calc_pe_brute(x, types, L, &potential, cutoff);
        
        int n_accept_translate = 0;
        int n_accept_swap = 0;
        int time_since=0;
        //make sure we clear out the visited histogram
        visited.clear_all();
        
        bool is_flat = false;

        int time_local = 0;
        while(is_flat == false)
        {

            if(check_max_time == true)
            {
                if(time_current> max_time)
                    is_flat = true;
            }
            /*if(time_local > 1000000)
            {
                std::cout << "*********************" << std::endl;
                std::cout << "things seem to be stuck, let's try stepping back" << std::endl;
                visited.clear_all();
                reset_configuration();
                ln_dos = ln_dos_old;
                time_local = 0;
            }*/
            int n_accept_local;
            
            //if(time_current == 155000)
            //    reset_configuration();
            
            if(translate_particles)
            {
                n_accept_local = translate_WL();
                n_accept_translate += n_accept_local;
                if(threads > 1)
                    reduce_dos();

            }
            nsq_neighbor_check_fast(x, nbr, cutoff, skin, L);

            //check to see if we need to rebuild the neighborlist
            if(swap_particles)
            {
                n_accept_local = swap_WL(); 
                n_accept_swap += n_accept_local;
                if(threads > 1)
                    reduce_dos();
            }
            //various outputs
            if(time_current%time_U_update == 0)
            {
                potential =0;
                calc_pe_brute(x, types, L, &potential, cutoff);
            }
            
          /*  if(time_current%time_lndos_swap == 0  && threads > 1)
            {
                reduce_dos();
            }
            else if(time_current%time_lndos_swap == 0  && threads == 1)
            {
                
                for(int i=0; i<visited.array.size(); i++)
                {
                    visited.array[i] = visited_local.array[i];
                    ln_dos.array[i] = ln_dos_local.array[i];
                }
            }*/

            
            if(time_current%time_thermo_output == 0)
            {
                double translate_probability = (double)n_accept_translate/(double)(time_thermo_output*x.size());
                double swap_probability = (double)n_accept_swap/(time_thermo_output*swap);
                
                print_thermo(thermoOut, time_current, rank, potential, N_particles, translate_probability, swap_probability);
                print_thermo(std::cout, time_current, rank, potential, N_particles, translate_probability, swap_probability);
                
                std::cout << "rank:\t" << rank;
                
                if(translate_particles)
                    std::cout << "\tdx:\t" << dx;
                if(swap_particles)
                    std::cout << "\tn swaps:\t" << swap;
                std::cout<< std::endl;
                
        
                loggerOut << time_current << "\t" << potential << std::endl;
                
                if(rank==0)
                {
                    char f1[1000];
                    sprintf(f1, "visited_temp.txt");
                    std::ofstream visited_f(f1);
                    visited.print(visited_f);
                    visited_f.close();
                    
                    sprintf(f1, "ln_dos_temp.txt");
                    std::ofstream ln_dos_f(f1);
                    ln_dos.print(ln_dos_f);
                    ln_dos_f.close();
                }
               /* if(time_current > time_thermo_output &&  translate_probability < 0.001 && time_since > 0)
                {
                    std::cout  <<"rank: " << rank << " reseting  configuration" << std::endl;

                    reset_configuration();
                    time_since = 0;
                }*/
                
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
                    }
                   /* if(swap_particles)
                    {
                        if(swap_probability < swap_target_prob && time_current > time_thermo_output)
                        {
                            swap++;
                            if(swap > swap_max)
                                swap = swap_max;
                        }
                        else if(swap_probability > swap_target_prob && time_current > time_thermo_output)
                        {
                            swap--;
                            if(swap < 1)
                            {
                                swap = 1;
                            }
                        }
                    }*/
                    n_accept_translate = 0;
                    n_accept_swap = 0;
                }
                
            }
            
            if(time_current%time_xyz_output == 0)
            {
                print_xyz(xyzOut, x, types, L);
            }
            
            if(time_current%time_check_flat == 0 && time_current > time_check_flat)
            {
                //if(threads > 1)
                //    check_reduce();
                double avg_flat=0, min_flat=0, max_flat=0, metric=0;
                bool check_flat =  check_flatness(visited, flatness, &avg_flat, &min_flat, &max_flat, &metric);
                //std::cout << "**** " << rank << "\t" << check_flat << std::endl;
                //if(rank ==0)
                //{
    
                    std::cout << std::setprecision(5) << "flatness min\t" << min_flat << "\tmax\t" << max_flat << "\tm\t" << metric << std::endl;
                    if(check_flat)
                    {
                        std::cout << "histogram converged with " << flatness << " for f = " << std::setprecision(32) << f << std::endl;
                        
                    }
                //}
                if(check_flat == true)
                {
                    is_flat = true;
                }
                
            }
            time_since++;
            time_current++;
            time_local++;
            
        }

        if(rank == 0)
        {
            char f1[1000];
            sprintf(f1, "visited_%.10f.txt", f);
            std::ofstream visitedOut(f1);
            visited.print(visitedOut);
            visitedOut.close();
            sprintf(f1, "ln_dos_%.10f.txt", f);
            std::ofstream dosOut(f1);
            ln_dos.print(dosOut);
            dosOut.close();
        }

        
    }
    
//    int translate_WL(coordlist_t& x, std::vector<int>& types, std::vector<neighbor> &nbr, double *L, double dx, double *potential_energy_total, histogram_1D& ln_dos, histogram_1D& visited_local, histogram_1D& ln_dos_local, double ln_factor, double min_energy, double max_energy)

    int translate_WL() 
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
                
                if(pe_new < U_max  && pe_new > U_min)
                {
                    /*if(pe_new < pe_min)
                        pe_min = pe_new;
                    */
                    double ln_dos_new = ln_dos.get_hist(pe_new);
                    double ln_dos_old = ln_dos.get_hist(pe_old);
                    
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
                        
                        visited.insert(pe_new);
                        ln_dos.insert(pe_new, ln_f);
                        
                        
                    }
                    else
                    {
                        visited.insert(pe_old);
                        ln_dos.insert(pe_old, ln_f);
                        
                    }
                }
                
            }
            
        }
        return n_accept;
        
        
    }
    int swap_WL()
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
        
        for(int i=0; i<swap; i++)
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
                        double pbc_new  = L[k]*anint(rij_new[k]*L_inv[k]);
                        
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
                
                if(pe_new < U_max  && pe_new > U_min)
                {
                    
                    /*if(pe_new < pe_min)
                        pe_min = pe_new;
                    */
                    double ln_dos_new = ln_dos.get_hist(pe_new);
                    double ln_dos_old = ln_dos.get_hist(pe_old);
                    
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
                        
                                                
                        visited.insert(pe_new);
                        ln_dos.insert(pe_new, ln_f);
                        
                        swapped[n1] = 1;
                        swapped[n2] = 1;
                    }
                    else
                    {
                        
                        visited.insert(pe_old);
                        ln_dos.insert(pe_old, ln_f);
                        
                    }
                }
            }
            
        }
        return n_accept;
        
        
    }


    
};


class WL_NVT_2D{
    
    
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
    double L[3];
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
    std::vector<int> types_init;
    
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
    
    std::ofstream loggerOut;
    
    bool setup_initialized;
    
    WL_NVT_2D()
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
        nn=0;
        time_current = 0;
    }
    
    void init_uniform_WL(double U_min_t, double U_max_t, double U_bin_t, double flatness_t, double M_min_t, double M_max_t, double M_bin_t)
    {
       /* U_min = U_min_t;
        U_max = U_max_t;
        U_bin = U_bin_t;*/
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
    void init_y(double M_min_t, double M_max_t, double M_bin_t)
    {
        M_min = M_min_t;
        M_max = M_max_t;
        M_bin = M_bin_t;
        visited.init_y(M_min, M_max, M_bin);
        ln_dos.init_y(M_min, M_max, M_bin);
        
    }
    
    void push_back_x(double U_min_t, double U_max_t, double U_bin_t)
    {
        visited.push_back_x(U_min_t, U_max_t, U_bin_t);
        ln_dos.push_back_x(U_min_t, U_max_t, U_bin_t);
    }
    
    void setup_histograms(double flatness_t)
    {
        visited.setup();
        ln_dos.setup();
        visited.check();
        ln_dos.check();
        
        flatness = flatness_t;
        WL_initialized = true;

    }
    
    void init_translate(double dx_t, double dx_min_t, double dx_max_t, double dx_target_prob_t)
    {
        dx = dx_t;
        dx_max = dx_max_t;
        dx_min = dx_min_t;
        dx_target_prob = dx_target_prob_t;
        
        translate_particles = true;
    }
    
    void init_system(double *L_t,  short unsigned int sseed_t)
    {
        for(int k=0; k<3; k++)
            L[k] = L_t[k];
        
        sseed = sseed_t;
        system_initialized = true;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        MPI_Comm_size(MPI_COMM_WORLD,&threads);
    }
    
    void init_swap(int swap_t, int swap_max_t, double swap_target_prob_t)
    {
        swap = swap_t;
        swap_max = swap_max_t;
        swap_target_prob = swap_target_prob_t;
        swap_particles = true;
        
        
    }
    void init_cswap(int swap_t, int swap_max_t, double swap_target_prob_t)
    {
        if(vol_change == false)
        {
            cswap = swap_t;
            cswap_max = swap_max_t;
            cswap_target_prob = swap_target_prob_t;
            cswap_particles = true;
        }
        else
        {
            std::cerr << "the 2nd WL dimension can either be a volume change or a composition change, not both at this time" << std::endl;
            assert(vol_change == false);
        }
        
    }
    
    void init_vol_change(int vol_t, int vol_max_t, double vol_target_prob_t, double ln_vol_change_max_t)
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
    
    
    void init_nlist(double cutoff_t, double skin_t)
    {
        cutoff = cutoff_t;
        skin = skin_t;
        
        nlist_initialized = true;
        
    }
    void set_configuration(coordlist_t x_t, std::vector<int> types_t)
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
                
                types.push_back(types_t[i]);
                types_init.push_back(types_t[i]);
                
            }
        }
        else
        {
            for(int i=0; i< x_t.size(); i++)
            {
                x[i] = x_t[i];
                x_init[i] = x_t[i];
                types[i] = types_t[i];
                types_init[i] = types_t[i];
            }
        }
        configuration_initialized = true;
    }
    void reset_configuration()
    {
        assert(x.size() == x_init.size() );
        
        for(int i=0; i<x.size(); i++)
        {
            x[i] = x_init[i];
            types[i] = types_init[i];
        }
        nsq_neighbor_rebuild(x, nbr, cutoff, skin, L);
        
        if(cswap_particles == true)
            count_type0();
        
    }
    
    
    
    void get_configuration(coordlist_t& x_t, std::vector<int>& types_t)
    {
        if(x.size() != x_t.size())
        {
            x_t.clear();
            types_t.clear();
            
            for(int i=0; i< x.size(); i++)
            {
                x_t.push_back(x[i]);
                types_t.push_back(types[i]);
            }
        }
        else
        {
            for(int i=0; i< x_t.size(); i++)
            {
                x_t[i] = x[i];
                types_t[i] = types[i];
            }
        }
    }
    void setup()
    {
        assert(dx < skin/2.0);
        for(int k=0; k<3; k++)
            assert(cutoff <= L[k]/2.0);
        assert(system_initialized);
        assert(nlist_initialized);
        assert(configuration_initialized);
        assert(WL_initialized);
        assert(x.size() > 0);
        nsq_neighbor_init(x, nbr, cutoff, skin, L);
        
        potential =0;
        calc_pe_brute(x, types, L, &potential, cutoff);
        if(translate_particles)
            nn+=x.size();
        
        char f1[1000];
        sprintf(f1, "pe_rank%d.txt", rank);
        loggerOut.open(f1);
        
        setup_initialized = true;
        
        
    }
    int get_threads()
    {
        MPI_Comm_size(MPI_COMM_WORLD,&threads);
        return threads;
    }
    int get_rank()
    {
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        return rank;
    }
    
    void init_ln_dos_from_file(const char* filename_ln_dos)
    {
        ln_dos.read_from_file(filename_ln_dos, get_threads());

       /* if(WL_initialized == true)
        {
            ln_dos.read_from_file2(filename_ln_dos, get_threads());
        }
        else
        {
            std::cerr<< "please define energy range and binsize first in init_WL" << std::endl;
            assert(WL_initialized==true);
        }*/
    }
    
    
    
    void reduce_dos()
    {
        
        
        visited.splat(0);
        ln_dos.splat(0);
        MPI_Allreduce(visited.array_local,visited.array,visited.size(),MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(ln_dos.array_local,ln_dos.array,ln_dos.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        //MPI::COMM_WORLD.Allreduce(visited.array_local,visited.array,visited.size(),MPI::DOUBLE,MPI::SUM);
        //MPI::COMM_WORLD.Allreduce(ln_dos.array_local,ln_dos.array,ln_dos.size(),MPI::DOUBLE,MPI::SUM);
        
        
        
    }
    void check_reduce()
    {
        double sum_S = 0;
        double sum_V = 0;
        
        for(int i=0; i<ln_dos.size(); i++)
        {
            sum_S+=ln_dos.array[i]/ln_dos.size();
            sum_V+=visited.array[i]/visited.size();
        }
        std::cout <<"check reduce:\t" << rank << "\t" << sum_S << "\t" << sum_V << std::endl;
    }
    
    void calc_U_min(double f_t, int run_time, int time_U_update, int time_check_flat, int time_thermo_output, bool adjust)
    {
        
        time_current = 0;
        assert(setup_initialized == true);
        
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
                if(threads > 1)
                    reduce_dos();
                
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
            
            if(time_current%time_thermo_output == 0)
            {
                std::cout << "M_current: " << M_current << std::endl;
                
                double translate_probability = (double)n_accept_translate/(double)(time_thermo_output*x.size());
                double cswap_probability = (double)n_accept_cswap/(time_thermo_output*cswap);
                
                print_thermo(std::cout, time_current, rank, potential, N_particles, translate_probability, cswap_probability);
                
                std::cout << "rank:\t" << rank;
                
                if(translate_particles)
                    std::cout << "\tdx:\t" << dx;
                if(cswap_particles)
                    std::cout << "\tn swaps:\t" << cswap;
                std::cout<< std::endl;
                
                
                loggerOut << time_current << "\t" << potential << std::endl;
                
                if(rank==0)
                {
                    char f1[1000];
                    sprintf(f1, "visited_temp_converge.txt");
                    std::ofstream visited_f(f1);
                    visited.print(visited_f);
                    visited_f.close();
                }
                
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
                }
                
            }
            
            
            if(time_current%time_check_flat == 0 && time_current > time_check_flat)
            {
                
                double avg_flat=0, min_flat=0, max_flat=0, metric=0;
                bool check_flat =  check_flatness(visited, flatness, &avg_flat, &min_flat, &max_flat, &metric);
                //std::cout << "**** " << rank << "\t" << check_flat << std::endl;
                if(rank ==0)
                {
                    
                    std::cout << std::setprecision(5) << "flatness min\t" << min_flat << "\tmax\t" << max_flat << "\tm\t" << metric << std::endl;
                    if(check_flat == true)
                    {
                        std::cout << "histogram converged with " << flatness << " for f = " << std::setprecision(32) << f << std::endl;
                        
                    }
                }
                if(check_flat == true)
                {
                    is_flat = true;
                }
                
            }
            time_since++;
            time_current++;
            time_local++;
            
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
        
    }
    
    void run(double f_t, int time_U_update, int time_check_flat, int time_xyz_output, std::ostream& xyzOut, int time_thermo_output, std::ostream& thermoOut, bool adjust)
    {
        
        time_current = 0;
        assert(setup_initialized == true);
        
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
            
            /*if(time_local > 1000000)
            {
                std::cout << "*********************" << std::endl;
                std::cout << "things seem to be stuck, let's try stepping back" << std::endl;
                visited.clear_all();
                reset_configuration();
                ln_dos = ln_dos_old;
                time_local = 0;
            }*/
            int n_accept_local;

            
            if(translate_particles)
            {
                n_accept_local = translate_WL();
                n_accept_translate += n_accept_local;
                if(threads > 1)
                    reduce_dos();
                
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
            
            if(time_current%time_thermo_output == 0)
            {
                std::cout << "M_current: " << M_current << std::endl;

                double translate_probability = (double)n_accept_translate/(double)(time_thermo_output*x.size());
                double cswap_probability = (double)n_accept_cswap/(time_thermo_output*cswap);
                
                print_thermo(thermoOut, time_current, rank, potential, N_particles, translate_probability, cswap_probability);
                print_thermo(std::cout, time_current, rank, potential, N_particles, translate_probability, cswap_probability);
                
                std::cout << "rank:\t" << rank;
                
                if(translate_particles)
                    std::cout << "\tdx:\t" << dx;
                if(cswap_particles)
                    std::cout << "\tn swaps:\t" << cswap;
                std::cout<< std::endl;
                
                
                loggerOut << time_current << "\t" << potential << std::endl;
                
                if(rank==0)
                {
                    char f1[1000];
                    sprintf(f1, "visited_temp.txt");
                    std::ofstream visited_f(f1);
                    visited.print(visited_f);
                    visited_f.close();
                    
                    sprintf(f1, "ln_dos_temp.txt");
                    std::ofstream ln_dos_f(f1);
                    ln_dos.print(ln_dos_f);
                    ln_dos_f.close();
                }
               /* if(time_current > time_thermo_output &&  translate_probability < 0.001 && time_since > 0)
                {
                    std::cout  <<"rank: " << rank << " reseting  configuration" << std::endl;
                    
                    reset_configuration();
                    time_since = 0;
                }*/
                
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
                }
                
            }
            
            if(time_current%time_xyz_output == 0)
            {
                print_xyz(xyzOut, x, types, L);
            }
            
            if(time_current%time_check_flat == 0 && time_current > time_check_flat)
            {
                //if(threads > 1)
                //    check_reduce();
                double avg_flat=0, min_flat=0, max_flat=0, metric=0;
                bool check_flat =  check_flatness(visited, flatness, &avg_flat, &min_flat, &max_flat, &metric);
                //std::cout << "**** " << rank << "\t" << check_flat << std::endl;
                //if(rank ==0)
                //{
                
                std::cout << std::setprecision(5) << "flatness min\t" << min_flat << "\tmax\t" << max_flat << "\tm\t" << metric << std::endl;
                if(check_flat)
                {
                    std::cout << "histogram converged with " << flatness << " for f = " << std::setprecision(32) << f << std::endl;
                    
                }
                //}
                if(check_flat == true)
                {
                    is_flat = true;
                }
                
            }
            time_since++;
            time_current++;
            time_local++;
            
        }
        
        if(rank == 0)
        {
            char f1[1000];
            sprintf(f1, "visited_%.10f.txt", f);
            std::ofstream visitedOut(f1);
            visited.print(visitedOut);
            visitedOut.close();
            sprintf(f1, "ln_dos_%.10f.txt", f);
            std::ofstream dosOut(f1);
            ln_dos.print(dosOut);
            dosOut.close();
        }
        
        
    }
    
    //    int translate_WL(coordlist_t& x, std::vector<int>& types, std::vector<neighbor> &nbr, double *L, double dx, double *potential_energy_total, histogram_1D& ln_dos, histogram_1D& visited_local, histogram_1D& ln_dos_local, double ln_factor, double min_energy, double max_energy)
    
    int translate_WL()
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
    int swap_WL()
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
        
        for(int i=0; i<swap; i++)
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
                        double pbc_new  = L[k]*anint(rij_new[k]*L_inv[k]);
                        
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
    int count_type0()
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
    int cswap_WL()
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
        
        for(int i=0; i<N_particles; i++)
        {
            int type0_temp = M_current;
            int n1,t1;
            
            n1 = i;
            t1 = types[n1];
      
        /*for(int i=0; i<cswap; i++)
        {
            bool picking = true;
            int type0_temp = M_current;
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
            */
           
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
        return n_accept;
        
        
    }
    double calc_vol()
    {
        M_current = L[0]*L[1]*L[2];
        return M_current;
    }
    
    
    
    int vol_WL()
    {
        return 0;
    }
    /*
        //this will calculate the volume and set M_current to this volume
        calc_vol();
        
        //precalculate 1/L
        
        double L_inv[3];
        for(int k=0; k<3; k++)
            L_inv[k] = 1.0/L[k];

        //vol corresponds to the number of attempted volume changes per cycle
        //in WL we could probably do a bunch of them each cycle
        //since ultimately we just want to sample different configurations
        for(int i=0; i<vol; i++)
        {
            double potential_energy_new = 0.0;
            double potential_energy_old = 0.0;
            bool overlap1 = false;
            bool overlap2 = false;
            
        }
        return 0;
    }
  /*
        //precalculate 1/L
        double L_inv[3];
        for(int k=0; k<3; k++)
            L_inv[k] = 1.0/L[k];
        
        calc_vol();
        
        
        double pot_temp = 0;
        
        int n_accept=0;
        
    
        for(int i=0; i<vol; i++)
        {            
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
                        double pbc_new  = L[k]*anint(rij_new[k]*L_inv[k]);
                        
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
        
        
    }*/

    
    
};
