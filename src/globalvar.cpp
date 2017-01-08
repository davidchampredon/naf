//
//  globalvar.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-16.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "globalvar.h"


// ==== Random seeds ====

uint	_RANDOM_SEED	= 12345;

std::mt19937_64	_RANDOM_GENERATOR(_RANDOM_SEED);
std::mt19937_64	_RANDOM_GENERATOR_INTERV(4321);

