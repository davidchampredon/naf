//
//  build_world.h
//  naf
//
//  Created by David CHAMPREDON on 2016-08-02.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#ifndef __naf__build_world__
#define __naf__build_world__

#include <stdio.h>

#include "socialPlace.h"
#include "dcTools.h"
#include "discrete_prob_dist.h"
#include "globalvar.h"
#include "individual.h"
#include "simulation.h"

using namespace std;

vector<individual> create_individuals(uint n);

vector<areaUnit> create_area_unit(const vector<ID>& id_au,
                                  const vector<string>& name_au,
                                  ID id_region,
                                  string regionName);


void keep_indiv_with_household(vector<individual>& x);

//vector<socialPlace> create_households_size(uint num_hh,
//                                           discrete_prob_dist<uint> size_distrib,
//                                           areaUnit AU,
//                                           vector<individual>& indiv);


vector<socialPlace> create_socialPlaces_size(SPtype sp_type,
                                             uint num_sp,
                                             discrete_prob_dist<uint> size_distrib,
                                             areaUnit AU,
                                             vector<individual>& indiv,
                                             float age_min = -999.999,
                                             float age_max = 999.99);

void assign_age_in_households(vector<socialPlace>& hh,
                              vector<individual>& indiv,
                              vector<vector<discrete_prob_dist<uint> > > age_distrib);


vector<socialPlace> build_world(vector<areaUnit> auvec,
                                discrete_prob_dist<uint> D_size_hh,     // Households sizes
                                vector< vector<discrete_prob_dist<uint> > > pr_age_hh,  // Age distribution inside households
                                discrete_prob_dist<uint> D_size_wrk,
                                discrete_prob_dist<uint> D_size_pubt,
                                discrete_prob_dist<uint> D_size_school,
                                vector<uint> n_hh,
                                vector<uint> n_wrk,
                                vector<uint> n_pubt,
                                vector<uint> n_school);


void move_all_indiv_to(SPtype sp_type,
                       world &w);

#endif /* defined(__naf__build_world__) */
