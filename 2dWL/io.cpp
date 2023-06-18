/*
 *  io.cpp
 *  simple_LJ
 *
 *  Created by Christopher Iacovella on 8/1/11.
 *  Copyright 2011 . All rights reserved.
 *
 */

#include "main.h"


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


std::ostream& operator << (std::ostream& os, coordlist_t& x)
{
	for (unsigned int i=0; i<x.size(); i++) {
		os << x[i] <<"\n";
	}
	return os;
}

std::ostream& operator << (std::ostream& os, coord_t& x)
{
	unsigned int sm1 = x.size()-1;
	for (unsigned int k=0; k<x.size(); k++) {
		os << x[k];
		if (k != sm1) {
			os <<"\t";
		}
	}
	return os;
}

std::istream& operator >> (std::istream& is, coordlist_t& x)
{
	return is;
}

std::istream& operator >> (std::istream& is, coord_t& x)
{
	return is;
}


void print(coordlist_t& x, std::ostream& os) 
{
	os << x << "\n";
}

void print(coord_t& x, std::ostream& os)
{
	os << x << "\n";
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

void print_xyz(const char* filename, coordlist_t& x, std::vector<int>& types, double *L)
{
	std::ofstream dataOut(filename);	
	
	print_xyz(dataOut, x, types, L);

}
void print_xyz(std::ostream& dataOut, coordlist_t& x, std::vector<int>& types, double *L)
{	
	dataOut << x.size() << std::endl;
	dataOut << "#\t" << L[0] << "\t" << L[1] << "\t" << L[2] << std::endl;
	for(int i=0; i<x.size(); i++)
	{
		dataOut << types[i] << "\t" << x[i] << std::endl;
	}
	
}

void load_xyz(const char* filename, coordlist_t& xyz, std::vector<int>& types)
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


void print_thermo(std::ostream& dataOut, int time, int rank, double potential, int N_particles, double translate_probability, double swap_probability)
{
    dataOut << time << "\trank\t" <<  rank << std::setprecision(5) <<"\tPE: " << potential << "\tPE/N: " << potential/(double)N_particles <<  std::setprecision(5) << "\ttrans. accept: " << translate_probability << "\tswap accept: " << swap_probability <<  std::endl;
}