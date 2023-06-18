%apply double *OUTPUT { double *average_flat, double *minimum_flat, double *maximum_flat}
%{
#include "flatness.hpp"
%}
%include "flatness.hpp"

