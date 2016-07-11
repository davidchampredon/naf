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
	
	
	
	
//	// Build world (test mode):
//	cout << "Generating world... "<<endl;
//	world W = test_world(sizereduction);
//	cout << "... world generated!"<<endl;
	
	
	// Build associated simulation:
	double horizon = 20;
	double sizereduction = 0.015; // Scale down world size compared to real world one
	
	Simulation sim;
	
	sim.build_test_world(sizereduction);
	
	auto t1 = std::chrono::system_clock::now();
	
	sim.set_horizon(horizon);
	
	if(sim.get_world().size()<2000) {
		sim.display_split_pop_linked();
		sim.display_split_pop_present();
	}
	
	// Seed infection(s) in world:
	
	vector<ID> id_pres = at_least_one_indiv_present(sim.get_world());
	
	
	vector<ID> sp_initially_infected {id_pres[0], id_pres[1]};
	vector<unsigned int> I0 {1,1};
	sim.seed_infection(sp_initially_infected, I0);
	
	if(sim.get_world().size()<2000) {
		sim.display_split_pop_present();
	}
	
	// Define model parameters:
	sim._modelParam.add_prm_double("proba_move", 0.90);
	sim._modelParam.add_prm_double("contact_rate", 3.0);
	
	// Run the simulation:
	sim.run();
	
	// timers:
	auto t2 = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = t2-t0;
	std::chrono::duration<double> elapsed_seconds2 = t2-t1;
	cout.precision(3);
	cout << "TOTAL TIME ELAPSED: "<< elapsed_seconds.count()/60.0 << " minutes" <<endl;
	cout << "Excluding pop generation: "<< elapsed_seconds2.count()/60.0 << " minutes" <<endl;
	
	return 0;
}
