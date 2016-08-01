//
//  tests.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-14.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include <stdio.h>
#include <chrono>

#include "tests.h"




Simulation test_naf(modelParam MP,
                    double horizon,
                    uint n_indiv,
                    uint i0,
                    const intervention &interv){
    
    // unpack parameters
    
    double dol_mean = MP.get_prm_double("dol_mean");
    double doi_mean = MP.get_prm_double("doi_mean");
    double doh_mean = MP.get_prm_double("doh_mean");
    bool debug_mode = MP.get_prm_bool("debug_mode");
    
    // Define the disease
    disease flu("Influenza", dol_mean, doi_mean, doh_mean);
    
    
    // Build simulation:
    
    Simulation sim;
    sim.set_modelParam(MP);
    sim.build_test_hospitalization(n_indiv);
    sim.set_horizon(horizon);
    sim.set_disease(flu);
    
    if(debug_mode){
        sim.display_split_pop_linked();
        sim.display_split_pop_present();
    }
    
    // Seed infection(s) in world:
    
    vector<ID> id_pres = at_least_one_indiv_present(sim.get_world());
    
    vector<ID> sp_initially_infected {id_pres[0]};
    vector<uint> I0 {i0};
    
    sim.define_all_id_tables();
    sim.seed_infection(sp_initially_infected, I0);
    
    sim.assign_hospital_to_individuals();
    
    sim.add_intervention(interv);
    
    // Run the simulation:
    sim.run();
    
    //sim.test();
    
    return sim;
}


void main_test_naf(){
    /// Main test
    
    double horizon = 300.0;
    modelParam MP;
    
    MP.add_prm_bool("debug_mode", true);
    
    MP.add_prm_string("dol_distrib", "exp");
    MP.add_prm_string("doi_distrib", "exp");
    MP.add_prm_string("doh_distrib", "exp");
    
    MP.add_prm_double ("dol_mean", 2.0);
    MP.add_prm_double ("doi_mean", 3.789);
    MP.add_prm_double ("doh_mean", 4.123);
    
    MP.add_prm_double ("proba_move", 1.0);
    
    MP.add_prm_double ("proba_death_prm_1", 0.59999999999999); // <-- TEST!
    MP.add_prm_double ("proba_death_prm_2", 0.85);
    MP.add_prm_double ("proba_death_prm_3", 0.80);
    
    MP.add_prm_bool   ("homogeneous_contact", false);
    MP.add_prm_double ("contact_rate", 2.0);
    MP.add_prm_uint   ("nt", 3);
    MP.add_prm_double ("asymptom_infectiousness_ratio", 0.8);
    MP.add_prm_double ("treat_doi_reduc", 1.123);
    MP.add_prm_double ("treat_reduc_infect_mean",0.1);
    
    MP.add_prm_double ("vax_imm_incr", 0.4);
    MP.add_prm_double ("vax_frail_incr",0.1);
    MP.add_prm_double ("vax_lag_full_efficacy", 99999);
    
    uint n_indiv = 500;
    uint i0 = 2;

    _RANDOM_GENERATOR.seed(123);
    
    // Define intervention
    float time_start = 2;
    float time_end = 999;
    float cvg_rate = 0.01;
    float cvg_max_prop = 0.30;
    intervention interv_test("vaccination",  // treatment  cure vaccination
                              "susceptible",  // symptomatic   susceptible
                              "interv_test",
                              time_start, time_end,
                              cvg_rate, cvg_max_prop);
    
    // Run simulations:
    
    Simulation sim1 = test_naf(MP, horizon, n_indiv, i0,interv_test);
    
    dcDataFrame pop_final = sim1.get_world()[2].export_dcDataFrame();
    pop_final.display();
}









Simulation test_SEIR_vs_ODE(modelParam MP,
							 double horizon,
							 uint n_indiv,
							 uint i0){
	
	// unpack parameters
	
    double dol_mean = MP.get_prm_double("dol_mean");
    double doi_mean = MP.get_prm_double("doi_mean");
    bool debug_mode = MP.get_prm_bool("debug_mode");
    
	// Define the disease
	disease flu("Influenza", dol_mean, doi_mean);
	
	// Build simulation:
	
	Simulation sim;
	sim.set_modelParam(MP);
	sim.build_single_world(n_indiv);
    sim.set_horizon(horizon);
	sim.set_disease(flu);
	
    if(debug_mode){
        sim.display_split_pop_linked();
        sim.display_split_pop_present();
    }
	
	// Seed infection(s) in world:
	
	vector<ID> id_pres = at_least_one_indiv_present(sim.get_world());
	
	vector<ID> sp_initially_infected {id_pres[0]};
	vector<uint> I0 {i0};
    
    sim.define_all_id_tables();
	sim.seed_infection(sp_initially_infected, I0);
    
	// Run the simulation:
	sim.run();
	
    sim.test();
    
	return sim;
}


void main_test_SEIR_vs_ODE(){
    /// function to be called in 'main.cpp'
    /// for testing SEIR vs ODE
    
    double horizon = 90.0;
    
    modelParam MP;
    
    MP.add_prm_bool("debug_mode", true);
    
    MP.add_prm_double("dol_mean", 2.0);
    MP.add_prm_double("doi_mean", 3.0);
    MP.add_prm_double("proba_move", 0.0);
    MP.add_prm_bool("homogeneous_contact", false);
    MP.add_prm_double("contact_rate", 1.0);
    MP.add_prm_uint("nt", 3);
    
    uint n_indiv = 1e4;
    uint i0 = 5;
    
    _RANDOM_GENERATOR.seed(123);
    Simulation sim1 = test_SEIR_vs_ODE(MP,horizon,n_indiv,i0);
    
    sim1.get_world()[0].export_dcDataFrame().display();
    sim1.timeseries().display();

}



Simulation test_move_2_sp(modelParam MP,
                          double horizon,
                          uint n_indiv,
                          uint i0){
    
    // unpack parameters
    
    double dol_mean = MP.get_prm_double("dol_mean");
    double doi_mean = MP.get_prm_double("doi_mean");
    bool debug_mode = MP.get_prm_bool("debug_mode");
    
    // Define the disease
    disease flu("Influenza", dol_mean, doi_mean);
    
    // Build simulation:
    
    Simulation sim;
    sim.set_modelParam(MP);

    sim.build_test_2_sp(n_indiv);
    
    sim.set_horizon(horizon);
    sim.set_disease(flu);
    
    if(debug_mode){
        sim.display_split_pop_linked();
        sim.display_split_pop_present();
    }
    
    // Seed infection(s) in world:
    
    vector<ID> id_pres = at_least_one_indiv_present(sim.get_world());
    
    vector<ID> sp_initially_infected {id_pres[0]};
    vector<uint> I0 {i0};
    
    sim.define_all_id_tables();
    sim.seed_infection(sp_initially_infected, I0);
    
    // Run the simulation:
    sim.run();
    
    //sim.test();
    
    return sim;
}


void main_test_move_2_sp(){
    /// Test moves of individuals between 2 social places
    
    
    double horizon = 30.0;
    
    modelParam MP;
    
    MP.add_prm_bool("debug_mode", true);
    
    MP.add_prm_double ("dol_mean", 2.0);
    MP.add_prm_double ("doi_mean", 3.0);
    MP.add_prm_double ("proba_move", 1.0);
    MP.add_prm_bool   ("homogeneous_contact", false);
    MP.add_prm_double ("contact_rate", 0.0);
    MP.add_prm_uint   ("nt", 3);

    uint n_indiv = 10;
    uint i0 = 2;
    
    
    _RANDOM_GENERATOR.seed(123);
    Simulation sim1 = test_move_2_sp(MP, horizon, n_indiv, i0);
}




Simulation test_hospitalization(modelParam MP,
                          double horizon,
                          uint n_indiv,
                          uint i0){
    
    // unpack parameters
    
    double dol_mean = MP.get_prm_double("dol_mean");
    double doi_mean = MP.get_prm_double("doi_mean");
    double doh_mean = MP.get_prm_double("doh_mean");
    bool debug_mode = MP.get_prm_bool("debug_mode");
    
    // Define the disease
    disease flu("Influenza", dol_mean, doi_mean, doh_mean);
    
    
    // Define intervention
    float time_start = 10;
    float time_end = 20;
    float cvg_rate = 0.05;
    float cvg_max_prop = 0.30;
    intervention interv_treat("cure",  // treatment  cure
                              "symptomatic",
                              "interv_treat",
                              time_start, time_end,
                              cvg_rate, cvg_max_prop);
    
    
    // Build simulation:
    
    Simulation sim;
    sim.set_modelParam(MP);
    sim.build_test_hospitalization(n_indiv);
    sim.set_horizon(horizon);
    sim.set_disease(flu);

    if(debug_mode){
        sim.display_split_pop_linked();
        sim.display_split_pop_present();
    }
    
    // Seed infection(s) in world:
    
    vector<ID> id_pres = at_least_one_indiv_present(sim.get_world());
    
    vector<ID> sp_initially_infected {id_pres[0]};
    vector<uint> I0 {i0};
    
    sim.define_all_id_tables();
    sim.seed_infection(sp_initially_infected, I0);
    
    sim.assign_hospital_to_individuals();
    
    sim.add_intervention(interv_treat);
    
    // Run the simulation:
    sim.run();
    
    //sim.test();
    
    return sim;
}


void main_test_hospitalization(){
    /// Test moves of individuals between 2 social places
    
    double horizon = 30.0;
    modelParam MP;
    
    MP.add_prm_bool("debug_mode", true);
    
    MP.add_prm_string("dol_distrib", "exp");
    MP.add_prm_string("doi_distrib", "exp");
    MP.add_prm_string("doh_distrib", "exp");
    
    MP.add_prm_double ("dol_mean", 2.0);
    MP.add_prm_double ("doi_mean", 3.789);
    MP.add_prm_double ("doh_mean", 4.123);

    MP.add_prm_double ("proba_move", 1.0);
    MP.add_prm_bool   ("homogeneous_contact", false);
    MP.add_prm_double ("contact_rate", 2.0);
    MP.add_prm_uint   ("nt", 3);
    MP.add_prm_double ("asymptom_infectiousness_ratio", 0.8);
    MP.add_prm_double ("doi_reduc_treat", 1.0);
    
    uint n_indiv = 200;
    uint i0 = 2;

    
    _RANDOM_GENERATOR.seed(123);
    Simulation sim1 = test_hospitalization(MP, horizon, n_indiv, i0);
    
    dcDataFrame pop_final = sim1.get_world()[2].export_dcDataFrame();
    pop_final.display();
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
	vector<uint> I0 {1,1};
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
	std::chrono::duration<double> elapsed_seconds = t2-t0;
	cout.precision(3);
	cout << "TIME 1: "<< elapsed_seconds.count()<< " secs" <<endl;
}


void test_random(){
	
	std::uniform_real_distribution<double> unif(0.0,1.0);
	for (int i=0; i<5; i++) {
		double y = unif(_RANDOM_GENERATOR);
		cout << "FCT-TEST-distrib = " << y <<endl;
	}
	
	
	
}











