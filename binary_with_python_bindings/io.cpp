/*
 *  io.cpp
 *  simple_LJ
 *
 *  Created by Christopher Iacovella on 8/1/11.
 *  Copyright 2011 . All rights reserved.
 *
 */

#include "WL.hpp"


bool check_line(std::string line)
{
    std::string space = " ";
    std::string tab = "\t";
    
    std::string comp1 = "";
    std::string comp2 = "\t";
    
    //a ghetto check to ensure that the line isn't just a bunch of white space or tabs
    for(int i=0; i<20; i++)
    {
        if(line == comp1)
        {
            return false;
        }
        comp1.append(space);
        if(line == comp2)
        {
            return false;
        }
        comp2.append(tab);
    }
    
    return true;
}



void load_raw(const char* filename, coordlist_t& x)
{
	std::ifstream file(filename);
	if (file.fail()) {
		std::cerr << "Error: load: can't open file " 
		<< filename << ".\n";
		exit(1);
	}
	
	std::string str;
	while (std::getline(file, str)) {
		std::istringstream iss(str.c_str());
		double temp;
		coord_t xi;
		while (iss >> temp) {
			xi.push_back(temp);
		}
		if (xi.size() != 0) {
			x.push_back(xi);
		}
	}
} 

void print_xyz_file(const char* filename, coordlist_t& x, types_t& types, coord_t L)
{
	std::ofstream dataOut(filename);	
	
	print_xyz(dataOut, x, types, L);

}
void print_xyz(std::ostream& dataOut, coordlist_t& x, types_t& types, coord_t L)
{	
	dataOut << x.size() << std::endl;
	dataOut << "#\t" << L[0] << "\t" << L[1] << "\t" << L[2] << std::endl;
	for(int i=0; i<x.size(); i++)
	{
		dataOut << types[i] << "\t" << x[i][0] << "\t" << x[i][1] << "\t" << x[i][2] << std::endl;
	}
	
}

void load_xyz(const char* filename, coordlist_t& xyz, types_t& types)
{
	int n_vectors=0;
	std::ifstream lp(filename);
	double x,y,z;
	int particle_type;
	
	std::cout << filename << std::endl;
	//make sure we can open the file
	assert(lp.good());
	
	xyz.clear();
	std::string temp_string, param;
	getline(lp, temp_string);
	std::istringstream iss(temp_string.c_str());
	iss >> n_vectors;
	
	getline(lp, temp_string);
	
	assert(n_vectors>0);
	
	int temp_type;
	
	for(int i=0; i<n_vectors; i++)
	{
		//make sure our file doesn't end early
		assert(lp.good());
		
		lp >> temp_type >> x >> y >> z;
        
		coord_t ctemp;
		ctemp.push_back(x);
		ctemp.push_back(y);
		ctemp.push_back(z);
		xyz.push_back(ctemp);
        types.push_back(temp_type);
	}
	
}


/*void print_thermo(std::ostream& dataOut, int time, int rank, double potential, int N_particles, double translate_probability, double swap_probability, double vol_probability)
{
    dataOut << time << "\trank\t" <<  rank << std::setprecision(5) <<"\tPE: " << potential << "\tPE/N: " << potential/(double)N_particles <<  std::setprecision(5) << "\ttrans. accept: " << translate_probability << "\tswap accept: " << swap_probability <<  "\tvol accept: " << vol_probability << std::endl;
}*/