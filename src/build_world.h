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
#include "simulator.h"

using namespace std;

vector<individual> create_individuals(uint n, uint first_id_indiv);

vector<areaUnit> create_area_unit(const vector<ID>& id_au,
                                  const vector<string>& name_au,
                                  ID id_region,
                                  string regionName);


void keep_indiv_with_household(vector<individual>& x);

void keep_indiv_with_relevant_links(vector<individual>& x);


vector<socialPlace> create_socialPlaces_size(SPtype sp_type,
                                             uint num_sp,
                                             uint first_id_sp,
                                             uint first_id_indiv,
                                             discrete_prob_dist<uint> size_distrib,
                                             areaUnit AU,
                                             vector<individual>& indiv,
                                             float age_min = -999.999,
                                             float age_max = 999.99);

vector<socialPlace> create_other_socialPlaces(uint num_sp,
                                              uint first_id_sp,
                                              discrete_prob_dist<uint> size_distrib,
                                              areaUnit AU);


void assign_age_in_households(vector<socialPlace>& hh,
                              vector<individual>& indiv,
                              vector<vector<discrete_prob_dist<uint> > > age_distrib);

/**
 * Build world before simulation runs.
 */
vector<socialPlace> build_world(vector<areaUnit> auvec,
                                vector<discrete_prob_dist<uint> > D_size_hh,     // Households sizes
                                vector< vector<discrete_prob_dist<uint> > > pr_age_hh,  // Age distribution inside households
                                
                                vector<discrete_prob_dist<uint> > D_size_wrk,
                                vector<discrete_prob_dist<uint> > D_size_pubt,
                                vector<discrete_prob_dist<uint> > D_size_school,
                                vector<discrete_prob_dist<uint> > D_size_hosp,
                                vector<discrete_prob_dist<uint> > D_size_other,
                                vector<uint> n_hh,
                                vector<uint> n_wrk,
                                vector<uint> n_pubt,
                                vector<uint> n_school,
                                vector<uint> n_hosp,
                                vector<uint> n_other);


/**
 * Assign schedule for all individual in the world.
 * For each individual, schedule chosen based on
 * age, employment status, etc.
 */
void assign_schedules(vector<socialPlace> & W,
                      const vector<schedule> & sched,
                      float prop_unemployed);

void assign_dox_distribution(vector<socialPlace> & W,
                             string dol_distrib,
                             string doi_distrib,
                             string doh_distrib);



void check_sp_integrity(const vector<socialPlace>& x);






#endif /* defined(__naf__build_world__) */
