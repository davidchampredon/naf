//
//  globalvar.h
//  naf
//
//  Created by David CHAMPREDON on 2016-07-16.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#ifndef __naf__globalvar__
#define __naf__globalvar__

#include <stdio.h>
#include <random>

// ==== Random seed ====

// The random number generator
// must be declared as a global variable
// in order to get the random seed
// initialization right for Rcpp.
// (I don't know why this has to work like that...)

extern std::mt19937_64	_RANDOM_GENERATOR;

// Need a secon, independent, random number generator
// in order not to interfer with simulation (e.g. mouvements, etc.)
// when individuals receive intervention before epidemic starts.
extern std::mt19937_64	_RANDOM_GENERATOR_INTERV;



// ==== File Directories ====


const std::string   _DIR_INPUT = "./INPUT/";


// ==== Global constants ====

const uint          _AGE_MAX = 99;

const float         AGE_YOUNG = 12.0;
const float         AGE_OLD   = 65.0;



#endif /* defined(__naf__globalvar__) */
