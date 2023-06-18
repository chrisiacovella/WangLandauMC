swig -python -c++ $1.i
#g++ -c -fPIC $1.cpp -o $1.o
g++ -c -fPIC *.cpp
g++ -c -fPIC $1_wrap.cxx -I$CPPFLAGS -L$LDFLAGS
ld -bundle -flat_namespace -undefined suppress -o _$1.so *.o
