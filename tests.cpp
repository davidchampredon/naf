//
//  tests.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-14.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include <stdio.h>




#include "tests.h"


void test_transmission(){
	
	auto t0 = std::chrono::system_clock::now();

	// Build associated simulation:
	
	double horizon = 50;
	Simulation sim;
	sim.build_single_world();
	sim.set_horizon(horizon);
	
	sim.display_split_pop_linked();
	sim.display_split_pop_present();
	
	auto t1 = std::chrono::system_clock::now();
	
	// Seed infection(s) in world:
	
	vector<ID> id_pres = at_least_one_indiv_present(sim.get_world());
	
	
	vector<ID> sp_initially_infected {id_pres[0], id_pres[1]};
	vector<unsigned int> I0 {1,1};
	sim.seed_infection(sp_initially_infected, I0);
	
	sim.display_split_pop_present();
	
	// Define model parameters:
	sim._modelParam.add_prm_double("proba_move", 0.90);
	sim._modelParam.add_prm_double("contact_rate", 1.5);
	
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

