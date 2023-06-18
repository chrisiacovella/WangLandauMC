/*
 *  main.h
 *  simple WL code
 *
 *  Created by Christopher Iacovella on 8/1/11.
 *  Copyright 2011. All rights reserved.
 *
 */

//#ifndef _WL_hpp
//#define _WL_hpp

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include <cmath>
#include <vector>
#include <fstream>

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <sstream>

#define kb 1

#define anint(x) ((x >= 0.5) ? (1.0) : (x <= -0.5) ? (-1.0) : (0.0))
#include "typedefs.hpp"
#include "average.hpp"
#include "system_config.hpp"
#include "initialize.hpp"
#include "random.hpp"
#include "pair_LJ.hpp"
#include "neighbor.hpp"
#include "energy.hpp"
#include "io.hpp"
#include "metropolis.hpp"
#include "histogram.hpp"
#include "flatness.hpp"
#include "wanglandau.hpp"

//#endif