
/* WL.i */

%module WL

%include "std_vector.i"


namespace std {
    %template(types_t) vector< int >;
    %template(coord_t)  vector < double >;
    %template(coordlist_t) vector < vector < double> >;
}
%include typemaps.i

%include typedefs.i
%include average.i
%include system_config.i
%include random.i
%include neighbor.i
%include initialize.i
%include pair_LJ.i
%include energy.i
%include io.i
%include metropolis.i
%include histogram.i
%include flatness.i
%include wanglandau.i
