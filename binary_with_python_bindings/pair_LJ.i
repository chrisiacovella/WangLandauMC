%apply double *OUTPUT { double *pe }
%{
#include "pair_LJ.hpp"
%}
%include "pair_LJ.hpp"

