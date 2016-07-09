//
//  main.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-04.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <limits.h>
#include <ctime>

#include <random>

#include "individual.h"
#include "socialPlace.h"
#include "simulation.h"
#include "utils.h"
#include "probaDistribution.h"

//#include "schedule.h"

using namespace std;


/*
 TO DO:
 - create large number of individuals, with age distribution
 - create large number of social places, associate individuals
 -
 
 */



int main(int argc, const char * argv[]) {
	
	auto t0 = std::chrono::system_clock::now();
	
	
	// ----------------------------------------------
//	probaDistrib<double> p(val,pr);
//	
//	vector<double> x = p.sample(50,1234);
//
//	displayVector(x);
//	
//	
//
//	
//	exit(99);
	// ----------------------------------------------
	
	string region1 = "Halton";
	ID id_region1 = 1;
	
	areaUnit A1(1, "Oakville", id_region1, region1);
	areaUnit A2(2, "Burlington", id_region1, region1);
	areaUnit A3(3, "Hamilton", id_region1, region1);
	
	vector<areaUnit> A {A1,A2,A3};
	
	// Schedule
	vector<double> timeslice {2.0/24, 8.0/24, 2.0/24, 12.0/24};
	vector<SPtype> worker_sed  {SP_pubTransp, SP_workplace, SP_pubTransp,SP_household};
	vector<SPtype> worker_trav {SP_pubTransp, SP_workplace, SP_other,SP_household};
	
	schedule sched_worker_sed(worker_sed, timeslice, "worker_sed");
	schedule sched_worker_trav(worker_trav, timeslice, "worker_trav");
	vector<schedule> sched {sched_worker_sed,sched_worker_trav};
	
	// test
	vector<individual> many_indiv = build_individuals(10000, sched);
	
	vector<SPtype> spt {SP_pubTransp, SP_workplace}; //, SP_other,SP_household, SP_school};
	vector<unsigned int> num_sp {50,80};
	
	probaDistrib<unsigned int> p_pubTransp({20,30,60},{0.6,0.3,0.1});
	probaDistrib<unsigned int> p_workPlace({5,20,50},{0.5,0.4,0.1});
	vector<probaDistrib<unsigned int> > p_size {p_pubTransp,p_workPlace};
	
	world Z = build_world_simple(spt, num_sp, p_size, many_indiv, A);
	
	Simulation toto(Z,20);
	toto.display_split_pop_linked();
	
	exit(99);
	
	double horiz = 20;
	unsigned int N_sp		= 7000;
	unsigned int N_indiv	= 50000;
	
	cout <<endl<<  " Generating population ..." << endl;
	world W = build_world_random(N_sp, A);
	populate_random_with_indiv(W, N_indiv, sched);
	cout << " done."<<endl;
	if(N_sp<1000) displayPopulationSize(W);
	
	auto t1 = std::chrono::system_clock::now();

	Simulation sim(W, horiz);
	
	sim.seed_infection({0,1}, {2,1});
	
	if(N_sp<1000) sim.display_split_pop_present();
	
	sim._modelParam.add_prm_double("proba_move", 0.90);
	sim._modelParam.add_prm_double("contact_rate", 3.0);

	sim.run();
	
	auto t2 = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = t2-t0;
	std::chrono::duration<double> elapsed_seconds2 = t2-t1;
	cout.precision(3);
	cout << "TOTAL TIME ELAPSED: "<< elapsed_seconds.count()/60.0 << " minutes" <<endl;
	cout << "Except pop generation: "<< elapsed_seconds2.count()/60.0 << " minutes" <<endl;
	return 0;
	}
