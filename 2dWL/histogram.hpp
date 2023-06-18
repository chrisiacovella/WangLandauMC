//
//  histogram.hpp
//
//
//  Created by Christopher Iacovella on 12/26/12.
//
//

#ifndef _histogram_hpp
#define _histogram_hpp

//basic class that will place data into a 1-D histogram

class histogram_1D{
public:
	//std::vector<double>array;
	double min;
	double max;
	double binsize;
    int bins;
	int total_entries;
    double *array;
    double *array_local;
    
    
    bool first_time;
	//sets the min, max and binsize
	histogram_1D(double mint, double maxt, double binsizet)
	{
		init(mint, maxt, binsizet);
	}
    histogram_1D()
    {
        //constructor that doesn't do anything
    }
    void init(double mint, double maxt, double binsizet)
    {
        min = mint; max = maxt;  binsize = binsizet;
        
        bins = (int)ceil((max-min)/binsize);
		
        total_entries =0;
        
        array = new double[bins];
        array_local = new double[bins];
        
        for(int i=0; i<bins; i++)
        {
            array[i] = 0.0;
            array_local[i] = 0.0;
        }
        
        first_time = true;
    }
    void splat(double value)
    {
        for(int i=0; i<bins; i++)
        {
            array[i] = value;
        }
    }
    
    void print(std::ostream& stream, double shift_bin)
	{
        double value = get_hist(shift_bin);
        
        double binsize2 = 0; //binsize/2.0;
		for(int i=0; i<bins; i++)
			stream << (i*binsize+min+binsize2) << "\t" << array[i]-value << std::endl;
	}
    
	void print(std::ostream& stream)
	{
        double binsize2 = 0; //binsize/2.0;
		for(int i=0; i<bins; i++)
			stream << (i*binsize+min+binsize2) << "\t" << std::setprecision(32) << array[i] << std::endl;
	}
	void print(const char* filename)
	{
        double binsize2 = 0; //binsize/2.0;
		std::ofstream stream(filename);
		for(int i=0; i<bins; i++)
			stream << (i*binsize+min+binsize2) << "\t" << std::setprecision(32) << array[i] << std::endl;
	}
	void print(const std::string& filename)
	{
        double binsize2 = 0; //binsize/2.0;
		std::ofstream stream(filename.c_str());
		for(int i=0; i<bins; i++)
			stream << (i*binsize+min+binsize2) << "\t" << std::setprecision(32) << array[i] << std::endl;
	}
    
    
    void read_from_file(const char* filename, int n_threads)
	{
		std::ifstream stream(filename);
        double temp_x, temp_y;
        //note this doesn't do any error checking!
		for(int i=0; i<bins; i++)
        {
			stream >> temp_x >> temp_y;
            array[i] = temp_y;
            //we'll just divide the ln_dos equally amongst the threads
            //that way the summation done during the mpi reduce is correct
            array_local[i] = temp_y/(double)n_threads;
            
        }
    }
    
    
	void insert(double pp)
	{
		double shifted = pp -min;
		int index = floor(shifted/binsize);
		if(index < bins)
		{
			array[index]++;
			array_local[index]++;
			total_entries++;
		}
	}
    void insert(double pp, double weight)
	{
		double shifted = pp -min;
		int index = floor(shifted/binsize);
		if(index < bins)
		{
			array[index]+=weight;
			array_local[index]+=weight;
			total_entries++;
		}
	}
    void scale(double pp, double factor)
	{
		double shifted = pp -min;
		int index = floor(shifted/binsize);
		if(index < bins)
		{
			array[index]*=factor;
		}
	}
	void normalize(double factor)
	{
		for(int i=0; i<bins; i++)
			array[i]/=factor;
	}
	void clear_global()
	{
		total_entries =0;
		for(int i=0; i<bins; i++)
			array[i]=0;
	}
    void clear_local()
	{
		total_entries =0;
		for(int i=0; i<bins; i++)
			array_local[i]=0;
	}
    void clear_all()
	{
		total_entries =0;
		for(int i=0; i<bins; i++)
        {
			array[i]=0;
   			array_local[i]=0;
        }
	}
    double get_hist(double pp)
    {
        double shifted = pp -min;
		int index = floor(shifted/binsize);
        if(index < bins)
            return array[index];
        else
            return 0;
    }
    
    double get_index(double pp)
    {
        double shifted = pp -min;
		int index = floor(shifted/binsize);
        
        return index;
    }
    double get_value(int index_t)
    {
        double val = index_t*binsize+min;
        return val;
    }
    
	//returns the total number of datapoints added to the histogram
	int size()
	{
		return bins;
	}
    
    int n_entries()
	{
		return total_entries;
	}

};

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
    
    histogram_2D()
    {
        initialized = false;
    }
    
    //let's set up a method that just makes a uniform array

    void init_uniform(double min_xt, double max_xt, double binsize_xt, double min_yt, double max_yt, double binsize_yt)
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
    void check()
    {
        assert(bins_y == container.size());
    }
    void init_y(double min_yt, double max_yt, double binsize_yt)
    {
        total_bins =0;
        //first, let's calculate the the number of "rows" we need
        min_y = min_yt; max_y = max_yt;  binsize_y = binsize_yt;
        bins_y = (int)ceil((max_y-min_y)/binsize_y);
        container.clear();
        initialized = false;

    }
    void push_back_x(double min_xt, double max_xt, double binsize_xt)
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
    
    void setup()
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
    std::vector<double> calculate_min()
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
    
    void insert(double pp_x, double pp_y, double weight)
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
    
    void insert(double pp_x, double pp_y)
    {
        insert(pp_x,  pp_y, 1.0);

    }
    
    void splat(double value)
    {
        for(int i=0; i<total_bins; i++)
        {
            array[i] = value;
        }
    }
    
    void print(std::ostream& stream)
	{
        
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
    void print(const char* filename)
    {
        std::ofstream stream(filename);
        print(stream);
    }
    void print(const std::string& filename)
    {
        std::ofstream stream(filename.c_str());
        print(stream);
    }

    //we should be able to parse this file without knowing any of the ranges
    //shouldn't be too hard
    /*
    void read_from_file2(const char* filename, int n_threads)
	{
        //first we are going to read the entire file in and determine the number of lines
        std::ifstream stream_check(filename);
        double temp_x, temp_y, temp_z;
        std::vector<int> ll;
        std::vector<double> low;
        std::vector<double> high;
        int total_length = 0;
        int total_y = 1;
        int total_x = 0;
        double old_y;

        double min_t;
        double max_t;
        container.clear();

        bool first_sweep = true;

        while(stream_check.good())
        {
            std::string line;
            while (std::getline(stream_check, line))
            {
                if(check_line(line))
                {
                    std::istringstream iss(line);
                    iss >> temp_x >> temp_y >> temp_z;
                    total_length++;
                    total_x++;
                    if(first_sweep == true)
                    {
                        old_y = temp_x;
                        first_sweep = false;
                        
                    }
                                       
                    if(temp_x != old_y)
                    {
                        old_y = temp_x;
                        ll.push_back(total_x);
                        total_x = 0;
                        total_y++;
                    }
                }
            }
        }
        std::cout << "total length " << total_length << "\ttotal rows: " << total_y << std::endl;
        
        
        array = new double[total_length];
        array_local = new double[total_length];
        
		std::ifstream stream(filename);
        while(stream.good())
        {
            stream >> temp_x >> temp_y >> temp_z;
            if(first_sweep == true)
            {
                old_y = temp_x;
                first_sweep = false;
            }
            if(temp_x == old_y)
            {
                container.
                
                
            }
            
            
            
        }
        double temp_x, temp_y, temp_z;
        //note this doesn't do any error checking!
        for(int j=0; j<bins_y; j++)
        {
            for(int i=0; i<container[j].bins; i++)
            {
                int ii = i+container[j].sum;
                stream >> temp_x >> temp_y >> temp_z;
                array[ii] = temp_z;
                //we'll just divide the ln_dos equally amongst the threads
                //that way the summation done during the mpi reduce is correct
                array_local[ii] = temp_z/(double)n_threads;
            }
        }
    }*/

    
    
    void read_from_file(const char* filename, int n_threads)
	{
		std::ifstream stream(filename);
        double temp_x, temp_y, temp_z;
        //note this doesn't do any error checking!
        for(int j=0; j<bins_y; j++)
        {
            for(int i=0; i<container[j].bins; i++)
            {
                int ii = i+container[j].sum;
                stream >> temp_x >> temp_y >> temp_z;
                array[ii] = temp_z;
                //we'll just divide the ln_dos equally amongst the threads
                //that way the summation done during the mpi reduce is correct
                array_local[ii] = temp_z/(double)n_threads;
            }
        }
    }
    void clear()
	{
		total_entries =0;
        for(int i=0; i<total_bins; i++)
			array[i]=0;
	}
    
    void clear_all()
	{
		total_entries =0;
        for(int i=0; i<total_bins; i++)
        {
			array[i]=0;
   			array_local[i]=0;
        }
	}
    void reset()
    {
        
        clear_all();
        container.clear();
        
    }
  
    double get_hist(double pp_x, double pp_y)
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

    int get_index(double pp_x, double pp_y)
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

    int size()
	{
		return total_bins;
	}
    
    int n_entries()
	{
		return total_entries;
	}
    double get_min(double pp_y)
    {
        double shifted_y = pp_y -min_y;
        int index_y = floor(shifted_y/binsize_y);

        return container[index_y].min;
    }
    double get_max(double pp_y)
    {
        double shifted_y = pp_y -min_y;
        int index_y = floor(shifted_y/binsize_y);
        
        return container[index_y].max;
    }
};

#endif
