//
//  tests.h
//  naf
//
//  Created by David CHAMPREDON on 2016-07-14.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#ifndef naf_tests_h
#define naf_tests_h

#include "individual.h"
#include "socialPlace.h"
#include "simulation.h"
#include "dcTools.h"
#include "probaDistribution.h"
#include "disease.h"
#include "globalvar.h"

Simulation test_transmission(modelParam MP,
							 double horizon,
							 unsigned int n_indiv,
							 unsigned int i0);

void test_move_transmission();


void test_rnd_eng();

void test_random();

#endif
