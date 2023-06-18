%apply double *OUTPUT { double *average_o, double *stdev_o }
%{
#include "average.hpp"
%}
%include "average.hpp"
