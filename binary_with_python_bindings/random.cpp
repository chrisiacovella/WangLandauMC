//
//  random.cpp
//  
//
//  Created by Chris Iacovella on 12/19/13.
//
//

#include "WL.hpp"

double rand_gaussian()
{
	double first,v1,v2,rsq,fac;
	static int save = 0;
	static double second;
	
	if (!save) {
		int again = 1;
		while (again) {
			v1 = 2.0*drand48()-1.0;
			v2 = 2.0*drand48()-1.0;
			rsq = v1*v1 + v2*v2;
			if (rsq < 1.0 && rsq != 0.0) again = 0;
                }
		fac = sqrt(-2.0*log(rsq)/rsq);
		second = v1*fac;
		first = v2*fac;
		save = 1;
	}
	else {
		first = second;
		save = 0;
	}
	return first;
}
//adding this because of strange issues importing the python random number generator
double rand_double()
{
    return drand48();
}