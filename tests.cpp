//
//  tests.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-14.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include <stdio.h>

#include "tests.h"


Simulation test_transmission(modelParam MP,
									   double horizon){
	
	auto t0 = std::chrono::system_clock::now();

	// unpack parameters
	
	double dol_mean			= MP.get_prm_double("dol_mean");
	double doi_mean			= MP.get_prm_double("doi_mean");
	unsigned int n_indiv	= MP.get_prm_uint("n_indiv");
	
	
	// Define the disease
	disease flu("Influenza", dol_mean, doi_mean);
	

	// Build simulation:
	
	Simulation sim;
	sim.set_modelParam(MP);
	sim.build_single_world(n_indiv);
	sim.set_horizon(horizon);
	sim.set_disease(flu);
	
	// Displays
	//sim.displayInfo_indiv();
	sim.display_split_pop_linked();
	sim.display_split_pop_present();
	
	auto t1 = std::chrono::system_clock::now();
	
	// Seed infection(s) in world:
	
	vector<ID> id_pres = at_least_one_indiv_present(sim.get_world());
	
	
	vector<ID> sp_initially_infected {id_pres[0], id_pres[1]};
	vector<unsigned int> I0 {1,1};
	sim.seed_infection(sp_initially_infected, I0);
	
	sim.display_split_pop_present();
	
	
	// Run the simulation:
	sim.run();
	
	// timers:
	auto t2 = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = t2-t0;
	std::chrono::duration<double> elapsed_seconds2 = t2-t1;
	cout.precision(3);
	cout << "TOTAL TIME ELAPSED: "<< elapsed_seconds.count()/60.0 << " minutes" <<endl;
	cout << "Excluding pop generation: "<< elapsed_seconds2.count()/60.0 << " minutes" <<endl;
	
	return sim;
}


void test_move_transmission(){
	
	auto t0 = std::chrono::system_clock::now();
	
	
	// Build associated simulation:
	double horizon = 20;
	double sizereduction = 0.002; // Scale down world size compared to real world one
	
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
}




void test_rnd_eng(){
	
	unsigned long N = 1e1;
	
	auto t0 = std::chrono::system_clock::now();
	
	auto t1 = std::chrono::system_clock::now();
	
	cout << endl;
	

	std::uniform_real_distribution<double> unif(0.0,1.0);
	
	for (int i=0; i<N; i++) {
		double y = unif(_RANDOM_GENERATOR);
		cout << "TEST-distrib = " << y <<endl;
	}
	cout << endl;
	


	std::uniform_real_distribution<double> unif2(0.0,1.0);
	
	for (int i=0; i<N; i++) {
		double y = unif2(_RANDOM_GENERATOR);
		cout << "TEST-distrib2 = " << y <<endl;
	}
	cout << endl;
	
	
	auto t2 = std::chrono::system_clock::now();
	
	
	
	std::chrono::duration<double> elapsed_seconds = t1-t0;
	std::chrono::duration<double> elapsed_seconds2 = t2-t1;
	cout.precision(3);
	cout << "TIME 1: "<< elapsed_seconds.count()<< " secs" <<endl;
	cout << "TIME 2: "<< elapsed_seconds2.count()<< " secs" <<endl;
}


void test_random(){
	
	std::uniform_real_distribution<double> unif(0.0,1.0);
	for (int i=0; i<5; i++) {
		double y = unif(_RANDOM_GENERATOR);
		cout << "FCT-TEST-distrib = " << y <<endl;
	}
	
	
	
}


