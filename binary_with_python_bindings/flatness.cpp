//
//  flatness.cpp
//  
//
//  Created by Chris Iacovella on 12/22/13.
//
//

#include "WL.hpp"

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
