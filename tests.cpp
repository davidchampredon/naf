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


void main_run_test(){
    
    // ================================================================
    //     MODEL PARAMETERS
    // ================================================================
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
    
    uint i0 = 2;
    
    _RANDOM_GENERATOR.seed(123);
    
    // Define intervention
    float time_start = 2;
    float time_end = 999;
    float cvg_rate = 0.01;
    float cvg_max_prop = 0.30;
    intervention interv("vaccination",  // treatment  cure vaccination
                        "susceptible",  // symptomatic   susceptible
                        "interv_test",
                        time_start, time_end,
                        cvg_rate, cvg_max_prop);

    // ================================================================
    //     WORLD PARAMETERS
    // ================================================================
    
    // Area Units
    vector<ID> id_au {0,1,}; //2,3};
    vector<string> name_au {"AUone", "AUtwo"}; //, "AUthree", "AUfour"};
    ID id_region = 0;
    string regionName = "RegionOne";
    vector<areaUnit> auvec = create_area_unit(id_au, name_au, id_region, regionName);
    
    // Households sizes
    vector<uint> hh_size {1,2,3};
    vector<double> hh_size_proba {0.3, 0.5, 0.2};
    discrete_prob_dist<uint> D_size_hh(hh_size, hh_size_proba);
    D_size_hh.display();
    
    // Age distribution inside households
    vector< vector<discrete_prob_dist<uint> > > pr_age_hh;
    
    discrete_prob_dist<uint> pr_age_hh_00({22,33,44},{0.2, 0.5, 0.3});
    discrete_prob_dist<uint> pr_age_hh_10({22,33,44},{0.2, 0.5, 0.3});
    discrete_prob_dist<uint> pr_age_hh_11({11,22,33,44},{0.2,0.1, 0.6, 0.1});
    discrete_prob_dist<uint> pr_age_hh_20({22,33,44},{0.2, 0.5, 0.3});
    discrete_prob_dist<uint> pr_age_hh_21({22,33,44},{0.2, 0.5, 0.3});
    discrete_prob_dist<uint> pr_age_hh_22({11,22,33,44},{0.8,0.1, 0.05, 0.05});
    
    vector<discrete_prob_dist<uint> > tmp;
    
    tmp.push_back(pr_age_hh_00);
    pr_age_hh.push_back(tmp);
    
    tmp.clear();
    tmp = {pr_age_hh_10,pr_age_hh_11};
    pr_age_hh.push_back(tmp);
    
    tmp.clear();
    tmp = {pr_age_hh_20,pr_age_hh_21,pr_age_hh_22};
    pr_age_hh.push_back(tmp);
    
    
    // Workplace sizes
    vector<uint> wrk_size {5,20,40,60};
    vector<double> wrk_size_proba {0.55, 0.3, 0.1, 0.05};
    discrete_prob_dist<uint> D_size_wrk(wrk_size, wrk_size_proba);
    D_size_wrk.display();
    
    // Public transport sizes
    vector<uint> pubt_size {20,50,60};
    vector<double> pubt_size_proba {0.3,0.4,0.3};
    discrete_prob_dist<uint> D_size_pubt(pubt_size, pubt_size_proba);
    D_size_pubt.display();
    
    // School sizes
    vector<uint> school_size {6, 8, 9};
    vector<double> school_size_proba {0.7,0.2,0.1};
    discrete_prob_dist<uint> D_size_school(school_size, school_size_proba);
    D_size_school.display();
    
    // Hospital sizes
    vector<uint> hosp_size {10000};
    vector<double> hosp_size_proba {1};
    discrete_prob_dist<uint> D_size_hosp(hosp_size, hosp_size_proba);
    D_size_hosp.display();
    
    // other public placed sizes
    vector<uint> other_size {10,50,100};
    vector<double> other_size_proba {0.4,0.4,0.2};
    discrete_prob_dist<uint> D_size_other(other_size, other_size_proba);
    D_size_other.display();
    
    
    // === Schedules definition ===
    
    // Define the time slices
    // must be same for all schedules and must sum up to 1.0
    vector<double> timeslice {1.0/24, 4.0/24, 4.0/24, 1.0/24, 2.0/24, 12.0/24};
    
    // type of schedules:
    vector<SPtype> worker_sed   {SP_pubTransp, SP_workplace, SP_workplace, SP_pubTransp, SP_other, SP_household};
    vector<SPtype> student      {SP_pubTransp, SP_school,    SP_school,    SP_pubTransp, SP_other, SP_household};
    vector<SPtype> unemployed   {SP_household, SP_other,     SP_other,     SP_other,     SP_other, SP_household};
    
    schedule sched_worker_sed(worker_sed, timeslice, "worker_sed");
    schedule sched_student(student, timeslice, "student");
    
    // Schedules used in the simulation:
    vector<schedule> sched {
        sched_worker_sed,
        sched_student
    };
    
    // number of each social place type
    vector<uint> n_hh       {50, 20};     //{5000,2500,1000,900};
    vector<uint> n_wrk      {30, 33};         //{10, 5, 3, 2};
    vector<uint> n_pubt     {40,30};          //{7,2,1,0};
    vector<uint> n_school   {20,20,1,1};
    vector<uint> n_hosp     {1,1,1,1};
    vector<uint> n_other    {100,50};       //{1000,500,250,200};
    
    
    // ================================================================
    //     RUN SIMULATION
    // ================================================================

    Simulation sim = run_test(auvec,
                              D_size_hh,
                              D_size_wrk,
                              D_size_pubt,
                              D_size_school,
                              D_size_hosp,
                              D_size_other,
                              pr_age_hh,
                              n_hh ,
                              n_wrk,
                              n_pubt ,
                              n_school,
                              n_hosp,
                              n_other,
                              sched ,
                              MP,
                              horizon,
                              i0,
                              interv);
    
    vector<dcDataFrame> df = export_dcDataFrame(sim.get_world());
    
}


Simulation run_test(vector<areaUnit> auvec,
                    discrete_prob_dist<uint> D_size_hh,
                    discrete_prob_dist<uint> D_size_wrk,
                    discrete_prob_dist<uint> D_size_pubt,
                    discrete_prob_dist<uint> D_size_school,
                    discrete_prob_dist<uint> D_size_hosp,
                    discrete_prob_dist<uint> D_size_other,
                    vector< vector<discrete_prob_dist<uint> > > pr_age_hh,
                    vector<uint> n_hh ,
                    vector<uint> n_wrk,
                    vector<uint> n_pubt ,
                    vector<uint> n_school,
                    vector<uint> n_hosp,
                    vector<uint> n_other,
                    vector<schedule> sched ,
                    modelParam MP,
                    double horizon,
                    uint i0,
                    const intervention &interv){
    
    
    
    Simulation sim;
    sim.set_modelParam(MP);
    
    sim.create_world(auvec,
                     D_size_hh, pr_age_hh,
                     D_size_wrk, D_size_pubt, D_size_school,D_size_school, D_size_other,
                     n_hh, n_wrk, n_pubt, n_school, n_hosp, n_other,
                     sched);
    
    sim.set_horizon(horizon);
    
    // Define the disease
    double dol_mean = MP.get_prm_double("dol_mean");
    double doi_mean = MP.get_prm_double("doi_mean");
    double doh_mean = MP.get_prm_double("doh_mean");
    bool debug_mode = MP.get_prm_bool("debug_mode");
    disease flu("Influenza", dol_mean, doi_mean, doh_mean);
    sim.set_disease(flu);
    
    if(debug_mode){
//        sim.display_split_pop_linked();
//        sim.display_split_pop_present();
    }
    
    // Seed infection(s) in world:
    vector<ID> id_pres = at_least_one_indiv_present(sim.get_world());
    vector<ID> sp_initially_infected {id_pres[0]};
    vector<uint> I0 {i0};
    sim.define_all_id_tables();
    sim.seed_infection(sp_initially_infected, I0);
    
    // Assign hospitals for all individuals:
    sim.assign_hospital_to_individuals();
    
    // Interventions:
    sim.add_intervention(interv);
    
    // Run the simulation:
    sim.run();
    
    return sim;
}




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
    
    uint n_indiv = 100;
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











