//
//  io.hpp
//  
//
//  Created by Chris Iacovella on 12/19/13.
//
//

#ifndef _io_hpp
#define _io_hpp

bool check_line(std::string line);
void load_raw(const char* filename, coordlist_t& x);
void print_xyz_file(const char* filename, coordlist_t& x, types_t& types, coord_t L);
void print_xyz(std::ostream& dataOut, coordlist_t& x, types_t& types, coord_t L);
void load_xyz(const char* filename, coordlist_t& xyz, types_t& types);

#endif
