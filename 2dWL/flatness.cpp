/*
 *  flatness.cpp
 *  simple_LJ
 *
 *  Created by Christopher Iacovella on 6/17/13.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "main.h"

void average_stdev(double *average_o, double *stdev_o, std::vector<double> data_list)
{
	double average=0.0, stdev=0.0;
	for(int i=0; i<data_list.size(); i++)
	{
		average+=data_list[i];
	}
	average = average/data_list.size();
	for(int i=0; i<data_list.size(); i++)
	{
		stdev+=(data_list[i]-average)*(data_list[i]-average);
	}
	stdev = sqrt(stdev/data_list.size());
	*average_o = average;
	*stdev_o = stdev;
}


//bool check_flatness(histogram_1D hist, double min, double max, double threshold, double *average_flat, double *minimum_flat, double *maximum_flat)
bool check_flatness(histogram_1D hist, double threshold, double *average_flat, double *minimum_flat, double *maximum_flat, double *metric)
{
    
    double diff = 1.0-threshold;
    double max_threshold = 1.0+diff;
    
    double average = 0;
    

    double count =0;
    double sum = 0;
    for(int i=0; i<hist.size(); i++)
    {
        sum += hist.array[i];
        count++;
    }

    average = sum/count;
    
    double avg, stdev;
   // average_stdev(&avg, &stdev, hist.array);
   // *metric = stdev/average;
    double variation = 0;
    if(average  != 0)
    {
        bool is_flat = true;
        double avg_flat = 0.0;
        double min_flat = 1.0;
        double max_flat = 0.0;
        int count = 0;
        
        int pad = 0;  //let's ignore the ends of the histogram when checking for flatness
        for(int i=pad; i<hist.size(); i++)
        {
            double fraction = hist.array[i]/average;
            if(fraction < min_flat)
                min_flat = fraction;
            if(fraction > max_flat)
                max_flat = fraction;
            
            if(fraction < threshold || fraction > max_threshold)
                is_flat = false;
            
            count++;
        }
        
        *metric = max_flat-1 - min_flat-1;
        *average_flat = average;
        *minimum_flat = min_flat-1;
        //*maximum_flat = max_flat;
        *maximum_flat = max_flat-1;
        return is_flat;
    }
    
    return false;
    
}
bool check_flatness(histogram_2D hist, double threshold, double *average_flat, double *minimum_flat, double *maximum_flat, double *metric)
{
    
    double diff = 1.0-threshold;
    double max_threshold = 1.0+diff;
    
    double average = 0;
    
    
    double count =0;
    double sum = 0;
    for(int i=0; i<hist.size(); i++)
    {
        sum += hist.array[i];
        count++;
    }
    
    average = sum/count;
    
    double avg, stdev;
    // average_stdev(&avg, &stdev, hist.array);
    // *metric = stdev/average;
    double variation = 0;
    if(average  != 0)
    {
        bool is_flat = true;
        double avg_flat = 0.0;
        double min_flat = 1.0;
        double max_flat = 0.0;
        int count = 0;
        
        int pad = 0;  //let's ignore the ends of the histogram when checking for flatness
        for(int i=pad; i<hist.size(); i++)
        {
            double fraction = hist.array[i]/average;
            if(fraction < min_flat)
                min_flat = fraction;
            if(fraction > max_flat)
                max_flat = fraction;
            
            if(fraction < threshold || fraction > max_threshold)
                is_flat = false;
            
            count++;
        }
        
        *metric = max_flat-1 - min_flat-1;
        *average_flat = average;
        *minimum_flat = min_flat-1;
        //*maximum_flat = max_flat;
        *maximum_flat = max_flat-1;
        return is_flat;
    }
    
    return false;
    
}


//
//bool check_flatness(histogram_2D hist, double threshold, double *average_flat, double *minimum_flat, double *maximum_flat)
//{
//    
//    double diff = 1.0-threshold;
//    double max_threshold = 1.0+diff;
//    
//    double average = 0;
//    
//    
//    double count =0;
//    double sum = 0;
//    for(int i=0; i<hist.hist.size(); i++)
//    {
//        for(int j=0; j<hist.hist[i].array.size(); j++)
//        {
//            sum += hist.hist[i].array[j];
//            count++;
//        }
//    }
//    
//    average = sum/count;
//    
//   // double avg, stdev;
//    //average_stdev(&avg, &stdev, hist.array);
//    
//    double variation = 0;
//    if(average  != 0)
//    {
//        bool is_flat = true;
//        double avg_flat = 0.0;
//        double min_flat = 1.0;
//        double max_flat = 0.0;
//        int count = 0;
//        for(int i=0; i<hist.hist.size(); i++)
//        {
//            for(int j=0; j<hist.hist[i].array.size(); j++)
//            {
//                
//                double fraction = hist.hist[i].array[j]/average;
//                if(fraction < min_flat)
//                    min_flat = fraction;
//                if(fraction > max_flat)
//                    max_flat = fraction;
//                
//                if(fraction < threshold || fraction > max_threshold)
//                    is_flat = false;
//            
//                count++;
//            }
//        }
//        
//        
//        *average_flat = average;
//        *minimum_flat = min_flat-1;
//        //*maximum_flat = max_flat;
//        *maximum_flat = max_flat-1;
//        return is_flat;
//    }
//    
//    return false;
//    
//}




