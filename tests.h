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


void main_test_naf();
Simulation test_naf(modelParam MP,
                    double horizon,
                    uint n_indiv,
                    uint i0,
                    const intervention &interv);



// --- One single world, SEIR model -------------------

Simulation test_SEIR_vs_ODE(modelParam MP,
                            double horizon,
                            uint n_indiv,
                            uint i0);

void main_test_SEIR_vs_ODE();
// ----------------------------------------------------



// --- Two SP, test movements -------------------
Simulation test_move_2_sp(modelParam MP,
                          double horizon,
                          uint n_indiv,
                          uint i0);

void main_test_move_2_sp();

// ----------------------------------------------------



// --- Test Hospitalization -------------------

Simulation test_hospitalization(modelParam MP,
                          double horizon,
                          uint n_indiv,
                          uint i0);
void main_test_hospitalization();
// ----------------------------------------------------




void test_move_transmission();
void test_rnd_eng();
void test_random();

#endif
