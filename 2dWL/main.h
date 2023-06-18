/*
 *  main.h
 *  simple_LJ
 *
 *  Created by Christopher Iacovella on 8/1/11.
 *  Copyright 2011. All rights reserved.
 *
 */
#include <mpi.h>

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include <cmath>
#include <vector>
#include <fstream>

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <sstream>

#define kb 1

#define anint(x) ((x >= 0.5) ? (1.0) : (x <= -0.5) ? (-1.0) : (0.0))

//define a custom vector type to make managing particle coordinates easier
typedef std::vector<double> coord_t;
typedef std::vector<coord_t> coordlist_t;


//define stream operators for the coord_t and coordlist class to make outputting data easier
std::ostream& operator << (std::ostream&, coordlist_t&);
std::ostream& operator << (std::ostream&, coord_t&);

std::istream& operator >> (std::istream&, coordlist_t&);
std::istream& operator >> (std::istream&, coord_t&);

void print(coordlist_t&, std::ostream&);
void print(coord_t&, std::ostream&);

void print_thermo(std::ostream& dataOut, int time, int rank, double potential, int N_particles, double translate_probability, double swap_probability);






//function I/O prototypes

void load_raw(const char* filename, coordlist_t& x);
void print_xyz(const char* filename, coordlist_t& x, std::vector<int>& types, double *L);
void print_xyz(std::ostream& dataOut, coordlist_t& x, std::vector<int>& types, double *L);
void init_system(coordlist_t& x, int N, double density, double *L, std::vector<int>& types);
void init_system_binary(coordlist_t& x, int N, int n_particles_0, double density, double *L, std::vector<int>& types);

void load_xyz(const char* filename, coordlist_t& xyz, std::vector<int>& types);

void check_cutoff(double cutoff, double skin, double *L);
bool check_line(std::string line);


#include "histogram.hpp"


//neighborlist class holds a particle id and the initial position of a particle
class neighbor{

public:
	std::vector <int> member;
	coordlist_t x_old;
    coord_t dx;
	
};


bool check_flatness(histogram_1D hist, double threshold, double *average_flat, double *minimum_flat, double *maximum_flat, double *metric);
bool check_flatness(histogram_2D hist, double threshold, double *average_flat, double *minimum_flat, double *maximum_flat, double *metric);



void nsq_neighbor_init(coordlist_t& x, std::vector<neighbor> &nbr, double cutoff, double skin, double *L);
void nsq_neighbor_check(coordlist_t& x, std::vector<neighbor> &nbr, double cutoff, double skin, double *L);
void nsq_neighbor_check_fast(coordlist_t& x, std::vector<neighbor> &nbr, double cutoff, double skin, double *L);
void nsq_neighbor_rebuild(coordlist_t& x, std::vector<neighbor> &nbr, double cutoff, double skin, double *L);

void calc_pe_brute(coordlist_t& x, std::vector<int>& types, double *L, double *potential_energy_total, double cutoff);

int pair_LJ_neighborlist(coordlist_t& x, std::vector<neighbor> &nbr, double sigma, double epsilon, double cutoff, double *L, double dx, double *potential_energy_total, double _beta, histogram_1D& vistied);

int pair_LJ_WL(coordlist_t& x, std::vector<int>& types, std::vector<neighbor> &nbr, double *L, double dx, double *potential_energy_total, histogram_1D& visited, histogram_1D& ln_dos, histogram_1D& visited_local, histogram_1D& ln_dos_local, double ln_factor, double min_energy, double max_energy, double min_energy_flat, double max_energy_flat);

/*int translate_NVT(coordlist_t& x, std::vector<int>& types, std::vector<neighbor> &nbr, double *L, double dx, double *potential_energy_total, double beta);
int swap_NVT(coordlist_t& x, std::vector<int>& types, std::vector<neighbor> &nbr, double *L, int n_swaps, double *potential_energy_total, double beta);


int translate_WL(coordlist_t& x, std::vector<int>& types, std::vector<neighbor> &nbr, double *L, double dx, double *potential_energy_total, histogram_1D& ln_dos, histogram_1D& visited_local, histogram_1D& ln_dos_local, double ln_factor, double min_energy, double max_energy);*/



int swap_LJ_WL(coordlist_t& x, std::vector<int>& types, std::vector<neighbor> &nbr, double *L, int n_swaps, double *potential_energy_total, histogram_1D& visited, histogram_1D& ln_dos, histogram_1D& visited_local, histogram_1D& ln_dos_local, double ln_factor, double min_energy, double max_energy, double min_energy_flat, double max_energy_flat);


#include "pair_LJ.hpp"
#include "metropolis.hpp"
#include "wanglandau.hpp"
