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

extern uint		_RANDOM_SEED;		// seed for random number generators
extern std::mt19937_64	_RANDOM_GENERATOR;




// ==== File Directories ====


const std::string    _DIR_INPUT = "./INPUT/";



#endif /* defined(__naf__globalvar__) */
