//
//  histogram.hpp
//
//
//  Created by Christopher Iacovella on 12/26/12.
//
//

#include "WL.hpp"


histogram_2D::histogram_2D()
{
    initialized = false;
}

    
void histogram_2D::init_uniform(double min_xt, double max_xt, double binsize_xt, double min_yt, double max_yt, double binsize_yt)
{
    
    if(initialized == false)
    {
        
        total_bins =0;
        //first, let's calculate the the number of "rows" we need
        min_y = min_yt; max_y = max_yt;  binsize_y = binsize_yt;
        bins_y = (int)ceil((max_y-min_y)/binsize_y);
        
        int bins_x = (int)ceil((max_xt-min_xt)/binsize_xt);
        
        container.clear();
        
        for(int i=0; i<bins_y; i++)
        {
            row_container temp;
            
            temp.min = min_xt;
            temp.max = max_xt;
            temp.binsize = binsize_xt;
            temp.bins = bins_x;
            temp.sum = 0;
            total_bins+=bins_x;
            temp.M = (double)i*binsize_y + min_y;
            
            if(i>0)
                temp.sum += container[i-1].bins + container[i-1].sum;
            
            container.push_back(temp);
        }
        
        total_entries =0;
        
        array = new double[total_bins];
        array_local = new double[total_bins];
        
        
        
        for(int i=0; i<total_bins; i++)
        {
            array[i] = 0.0;
            array_local[i] = 0.0;
        }
        
        first_time = true;
        initialized = true;
    }
    else
    {
        std::cerr << "histogram has already been initialized!\n" << std::endl;
        assert(initialized == false);
    }
}

void histogram_2D::read_histogram_params(const char* filename)
{
    std::ifstream stream(filename);
    
    //read back in the histogram parameter file
    
    if(initialized == false)
    {
        std::string temp_string;
        getline(stream, temp_string);
        
        std::istringstream iss(temp_string.c_str());
        int temp_value;
        iss >> temp_value;
        total_bins = temp_value;
        
        double min_yt, max_yt, binsize_yt;
        int bins_yt;
        
        std::istringstream iss2(temp_string.c_str());
        iss2 >> min_yt, max_yt, binsize_yt, bins_yt;
        
        min_y = min_yt;
        max_y = max_yt;
        binsize_y = binsize_yt;
        bins_y = bins_yt;
        
        
        
        container.clear();
        
        double min_xt, max_xt, binsize_xt;
        int bins_xt;
        
        for(int i=0; i<bins_y; i++)
        {
            row_container temp;
            
            std::istringstream iss3(temp_string.c_str());
            iss3 >> min_xt >> max_xt >> binsize_xt >> bins_xt;
            
            temp.min = min_xt;
            temp.max = max_xt;
            temp.binsize = binsize_xt;
            temp.bins = bins_xt;
            temp.sum = 0;
            temp.M = (double)i*binsize_y + min_y;
            
            if(i>0)
                temp.sum += container[i-1].bins + container[i-1].sum;
            
            container.push_back(temp);
            
        }
        
        
        
        
        total_entries =0;
        
        array = new double[total_bins];
        array_local = new double[total_bins];
        
        
        
        for(int i=0; i<total_bins; i++)
        {
            array[i] = 0.0;
            array_local[i] = 0.0;
        }
        
        first_time = true;
        initialized = true;
    }
    else
    {
        std::cerr << "histogram has already been initialized!\n" << std::endl;
        assert(initialized == false);
    }
}


void histogram_2D::check()
{
    assert(bins_y == container.size());
}
void histogram_2D::init_y(double min_yt, double max_yt, double binsize_yt)
{
    total_bins =0;
    //first, let's calculate the the number of "rows" we need
    min_y = min_yt; max_y = max_yt;  binsize_y = binsize_yt;
    bins_y = (int)ceil((max_y-min_y)/binsize_y);
    container.clear();
    initialized = false;
    
}
void histogram_2D::push_back_x(double min_xt, double max_xt, double binsize_xt)
{
    if(initialized == false)
    {
        
        
        int bins_x = (int)ceil((max_xt-min_xt)/binsize_xt);
        
        row_container temp;
        
        temp.M = (double)container.size()*binsize_y + min_y;
        temp.min = min_xt;
        temp.max = max_xt;
        temp.binsize = binsize_xt;
        temp.bins = bins_x;
        temp.sum = 0;
        total_bins+=bins_x;
        if(container.size() > 0)
            temp.sum += container[container.size()-1].bins + container[container.size()-1].sum;
        
        container.push_back(temp);
    }
    else
    {
        std::cerr << "histogram has already been initialized!\n" << std::endl;
        assert(initialized == false);
    }
}

void histogram_2D::setup()
{
    if(initialized==false)
    {
        
        array = new double[total_bins];
        array_local = new double[total_bins];
        
        for(int i=0; i<total_bins; i++)
        {
            array[i] = 0.0;
            array_local[i] = 0.0;
        }
        
        
        
        
        initialized = true;
    }
    else
    {
        std::cerr << "histogram has already been initialized!\n" << std::endl;
        assert(initialized == false);
    }
}
std::vector<double> histogram_2D::calculate_min()
{
    int count = 0;
    std::vector<double> temp_Umin;
    for(int j=0; j<bins_y; j++)
    {
        double lowest_U=container[j].max;
        for(int i=0; i<container[j].bins; i++)
        {
            if(array[count] > 0)
                if((i*container[j].binsize+container[j].min) < lowest_U)
                {
                    lowest_U = ((i)*container[j].binsize+container[j].min);
                }
            count++;
        }
        temp_Umin.push_back(lowest_U+5.0*container[j].binsize);
    }
    return temp_Umin;
}

void histogram_2D::insert(double pp_x, double pp_y, double weight)
{
    if(initialized == true)
    {
        //we first need to figure out which row we are in, index_y
        
        double shifted_y = pp_y -min_y;
        int index_y = floor(shifted_y/binsize_y);
        
        
        double shifted_x = pp_x -container[index_y].min;
        
        int index_x = floor(shifted_x/container[index_y].binsize);
        
        
        if(index_x < container[index_y].bins && index_y < bins_y)
        {
            int ii = index_x + container[index_y].sum;
            array[ii]+=weight;
            array_local[ii]+=weight;
            total_entries++;
        }
    }
    else
    {
        std::cerr << "histogram needs to be initialized!\n" << std::endl;
        assert(initialized == true);
        
    }
}

void histogram_2D::insert(double pp_x, double pp_y)
{
    insert(pp_x,  pp_y, 1.0);
    
}

void histogram_2D::splat(double value)
{
    for(int i=0; i<total_bins; i++)
    {
        array[i] = value;
    }
}

void histogram_2D::print_histogram(const char* filename)
{
    std::ofstream stream(filename);
    
    int count = 0;
    for(int j=0; j<bins_y; j++)
    {
        for(int i=0; i<container[j].bins; i++)
        {
            stream << (j*binsize_y+min_y)  << "\t" << (i*container[j].binsize+container[j].min) << "\t" << std::setprecision(32) <<  array[count] << std::endl;
            count++;
        }
    }
}


void histogram_2D::print_params(const char* filename)
{
    std::ofstream stream(filename);
    //we're going to output all the parameters we need to easily read the script back in
    
    stream << total_bins << "\t#total bins" << std::endl;
    stream << min_y << "\t" << max_y << "\t" << binsize_y << "\t" << bins_y << "\t#min_y max_y binsize_y bins_y" << std::endl;
    
    for(int i=0; i<container.size(); i++)
    {
        stream << container[i].min << "\t" << container[i].max << "\t" << container[i].binsize << "\t" << container[i].bins  << "\t#min_x max_x binsize_x bins_x" << std::endl;
    }
}



void histogram_2D::read_from_file(const char* filename, int n_threads)
{
    std::ifstream stream(filename);
    double temp_x, temp_y, temp_z;
    //note this doesn't do any error checking!
    for(int j=0; j<total_bins; j++)
    {
        stream >> temp_x >> temp_y >> temp_z;
        array[j] = temp_z;
        array_local[j] = temp_z/(double)n_threads;
        
    }
}
void histogram_2D::clear()
{
    total_entries =0;
    for(int i=0; i<total_bins; i++)
        array[i]=0;
}

void histogram_2D::clear_all()
{
    total_entries =0;
    for(int i=0; i<total_bins; i++)
    {
        array[i]=0;
        array_local[i]=0;
    }
}
void histogram_2D::reset()
{
    
    clear_all();
    container.clear();
    delete [] array;
    delete [] array_local;
}

double histogram_2D::get_hist(double pp_x, double pp_y)
{
    
    double shifted_y = pp_y -min_y;
    int index_y = floor(shifted_y/binsize_y);
    
    
    double shifted_x = pp_x -container[index_y].min;
    
    int index_x = floor(shifted_x/container[index_y].binsize);
    
    
    
    if(index_x < container[index_y].bins && index_y < bins_y)
    {
        int ii = index_x + container[index_y].sum;
        return array[ii];
    }
    else
    {
        return 0.0;
        
    }
}

int histogram_2D::get_index(double pp_x, double pp_y)
{
    double shifted_y = pp_y -min_y;
    int index_y = floor(shifted_y/binsize_y);
    
    
    double shifted_x = pp_x -container[index_y].min;
    
    int index_x = floor(shifted_x/container[index_y].binsize);
    
    
    
    if(index_x < container[index_y].bins && index_y < bins_y)
    {
        int ii = index_x + container[index_y].sum;
        return ii;
    }
    else
        return -1;
}
int histogram_2D::size()
{
    return total_bins;
}

int histogram_2D::n_entries()
{
    return total_entries;
}
double histogram_2D::get_min(double pp_y)
{
    double shifted_y = pp_y -min_y;
    int index_y = floor(shifted_y/binsize_y);
    
    return container[index_y].min;
}
double histogram_2D::get_max(double pp_y)
{
    double shifted_y = pp_y -min_y;
    int index_y = floor(shifted_y/binsize_y);
    
    return container[index_y].max;
}
