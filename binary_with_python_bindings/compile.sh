swig -python -c++ WL.i
#g++ -c -fPIC $1.cpp -o $1.o
g++ -c -fPIC -O2 *.cpp
g++ -c -fPIC WL_wrap.cxx -I/Users/cri/miniforge3/envs/WL/include/python3.11/ -L/Users/cri/miniforge3/envs/WL/lib/python3.11/
ld -bundle -flat_namespace -undefined suppress -o _$1.so *.o
