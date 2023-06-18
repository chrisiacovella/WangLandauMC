//
//  average.cpp
//  
//
//  Created by Chris Iacovella on 12/22/13.
//
//

#include "WL.hpp"


void average_stdev(double *average_o, double *stdev_o, coord_t data_list)
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