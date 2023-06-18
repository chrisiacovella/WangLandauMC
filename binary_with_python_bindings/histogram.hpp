//
//  histogram.hpp
//
//
//  Created by Christopher Iacovella on 12/26/12.
//
//

#ifndef _histogram_hpp
#define _histogram_hpp

//row container will hold all the information about a given row, but no data.
//... in this case, min and max energy, plus the size the array needs to be.
class row_container{
public:
    double min;
    double max;
    int bins;
    double binsize;
    double M;
    int sum;
};

//we need to have a 2d histogram that will allow a unique energy range for each of bins
//in the second dimension
//the easiest way to do this might be to make an array of 1d histograms
//but makes problems for MPI, since it's easiest to just use MPI_ALLREDUCE of a double array
//I think I just need to make up a vector that keeps track of what each bin in the array corresponds to


class histogram_2D{
public:
    std::vector<row_container> container;
    
    double min_y;
	double max_y;
	double binsize_y;
    
    int bins_y;
    
    int total_bins;
	int total_entries;

    double *array;
    double *array_local;
    
    bool first_time;
    bool initialized;
    
    histogram_2D();
    
    //let's set up a method that just makes a uniform array

    void init_uniform(double min_xt, double max_xt, double binsize_xt, double min_yt, double max_yt, double binsize_yt);
    void read_histogram_params(const char* filename);
    void check();
    void init_y(double min_yt, double max_yt, double binsize_yt);
    void push_back_x(double min_xt, double max_xt, double binsize_xt);
    void setup();
    std::vector<double> calculate_min();
    void insert(double pp_x, double pp_y, double weight);
    void insert(double pp_x, double pp_y);
    void splat(double value);
    void print_histogram(const char* filename);
    
    void print_params(const char* filename);
  

    void read_from_file(const char* filename, int n_threads);
    void clear();
    void clear_all();
    void reset();
    double get_hist(double pp_x, double pp_y);
    int get_index(double pp_x, double pp_y);
    int size();
    int n_entries();
    double get_min(double pp_y);
    double get_max(double pp_y);
};

#endif
