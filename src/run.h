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
#include "simulator.h"
#include "dcTools.h"
#include "discrete_prob_dist.h"
#include "disease.h"
#include "globalvar.h"
#include "build_world.h"
#include "dcMatrix.h"


Simulator run_stochWorld(vector<areaUnit> auvec,
                   vector<discrete_prob_dist<uint> > D_size_hh,
                   vector<discrete_prob_dist<uint> > D_size_wrk,
                   vector<discrete_prob_dist<uint> > D_size_pubt,
                   vector<discrete_prob_dist<uint> > D_size_school,
                   vector<discrete_prob_dist<uint> > D_size_hosp,
                   vector<discrete_prob_dist<uint> > D_size_other,
                   vector< vector<discrete_prob_dist<uint> > > pr_age_hh,
                   vector<uint> n_hh ,
                   vector<uint> n_wrk,
                   vector<uint> n_pubt ,
                   vector<uint> n_school,
                   vector<uint> n_hosp,
                   vector<uint> n_other,
                   float unemployed_prop,
                   vector<schedule> sched ,
                   modelParam MP,
                   double start_time,
                   double horizon,
                   uint i0,
                   const vector<intervention> &interv,
                   bool build_world_only);


Simulator run_detWorld(vector<areaUnit> auvec,
                       vector<vector<uint> > size_hh,
                       vector<vector<uint> > size_wrk,
                       vector<vector<uint> > size_pubt,
                       vector<vector<uint> > size_school,
                       vector<vector<uint> > size_hosp,
                       vector<vector<uint> > size_other,
                       vector< vector<discrete_prob_dist<uint> > > pr_age_hh,
                       float unemployed_prop,
                       vector<schedule> sched ,
                       modelParam MP,
                       double start_time,
                       double horizon,
                       uint i0,
                       const vector<intervention> &interv,
                       bool build_world_only);


void main_run();


// -----------------------------------------
void main_test_naf();
Simulator test_naf(modelParam MP,
                    double horizon,
                    uint n_indiv,
                    uint i0,
                    const intervention &interv);
// -----------------------------------------



// --- One single world, SEIR model -------------------

Simulator test_SEIR_vs_ODE(modelParam MP,
                            double horizon,
                            uint n_indiv,
                            uint i0);

void main_test_SEIR_vs_ODE();
// ----------------------------------------------------



// --- Two SP, test movements -------------------
Simulator test_move_2_sp(modelParam MP,
                          double horizon,
                          uint n_indiv,
                          uint i0);

void main_test_move_2_sp();

// ----------------------------------------------------



// --- Test Hospitalization -------------------

Simulator test_hospitalization(modelParam MP,
                          double horizon,
                          uint n_indiv,
                          uint i0);
void main_test_hospitalization();
// ----------------------------------------------------



void test_move_transmission();
void test_rnd_eng();
void test_random();


#endif
