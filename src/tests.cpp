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


void main_run(){
    
    // ================================================================
    //     MODEL PARAMETERS
    // ================================================================
    
    double start_time = -5.0;
    double horizon    = 300.0;
    
    modelParam MP;
    
    MP.add_prm_bool("debug_mode", false);
    
    MP.add_prm_string("dol_distrib", "lognorm");
    MP.add_prm_string("doi_distrib", "lognorm");
    MP.add_prm_string("doh_distrib", "exp");
    
    MP.add_prm_double ("dol_mean", 2.0);
    MP.add_prm_double ("doi_mean", 4.789);
    MP.add_prm_double ("doh_mean", 4.123);
    
    MP.add_prm_double ("dol_var", 2.1);
    MP.add_prm_double ("doi_var", 2.5);
    MP.add_prm_double ("doh_var", 0.2);
    
    MP.add_prm_double ("proba_move", 0.9);
    MP.add_prm_double ("proba_move_reduc_sympt", 0.15);

    MP.add_prm_double ("proba_hosp", 0.01);
    
    MP.add_prm_double ("proba_change_sp_other", 0.01);
    
    MP.add_prm_double ("proba_death_prm_1", 0.59999999999999); // <-- TEST!
    MP.add_prm_double ("proba_death_prm_2", 0.85);
    MP.add_prm_double ("proba_death_prm_3", 0.80);
    
    MP.add_prm_double("frailty_0", 0.60);
    MP.add_prm_double("frailty_min", 0.15);
    MP.add_prm_double("frailty_agemin", 30.0);
    MP.add_prm_double("frailty_agepivot", 60);
    MP.add_prm_double("frailty_pivot", 0.50);
    MP.add_prm_double("frailty_powerChild", 3);
    MP.add_prm_double("frailty_sd", 0.1);
    
    MP.add_prm_double("imm_hum_baseline", 0.1);
    MP.add_prm_double("imm_hum_agezero", 100);
    MP.add_prm_double("imm_hum_p", 2.0);
    
    MP.add_prm_double("imm_cell_max", 0.7);
    MP.add_prm_double("imm_cell_slope", 2);
    MP.add_prm_double("imm_cell_pivot", 20);
    
    
    MP.add_prm_bool   ("homogeneous_contact", false);
    MP.add_prm_double ("contact_rate_mean", 4.98);
    MP.add_prm_double ("contact_rate_stddev", 1);
    MP.add_prm_double ("contact_ratio_age_1_10", 2.0);
    MP.add_prm_double ("contact_ratio_age_10_16", 1.5);
    MP.add_prm_double ("contact_ratio_age_over_65", 0.8);
    MP.add_prm_double ("contact_ratio_sp_household", 2.2);
    MP.add_prm_double ("contact_ratio_sp_pubTransport", 1.75);
    
    MP.add_prm_double("contactAssort_lambda", 1.0/1.0);
    
    MP.add_prm_uint   ("nt", 3);
    
    MP.add_prm_double ("asymptom_infectiousness_ratio", 0.8);
    MP.add_prm_double ("treat_doi_reduc", 1.123);
    MP.add_prm_double ("treat_reduc_infect_mean",0.1);
    
    MP.add_prm_double ("vax_imm_hum_incr", 0.2);
    MP.add_prm_double ("vax_imm_cell_incr", 0.4);
    MP.add_prm_double ("vax_frail_incr",0.1);
    MP.add_prm_double ("vax_lag_full_efficacy", 99999);
    
    uint i0 = 10;
    
    _RANDOM_GENERATOR.seed(123);
    
    // Define intervention
    vector<intervention> interv_vec;
    
    float time_start   = -2;
    float time_end     = 999;
    float cvg_rate     = 0.02;
    float cvg_max_prop = 0.30;
    intervention interv1("vaccination",  // treatment  cure vaccination
                         "young_old",  // symptomatic   susceptible   young_old
                         "vax",
                         time_start, time_end,
                         cvg_rate, cvg_max_prop);
    
    float time_start2    = 0;
    float time_end2      = 999;
    float cvg_rate2      = 0.60;
    float cvg_max_prop2  = 0.99;
    intervention interv2("treatment",  // treatment  cure vaccination
                         "symptomatic",  // symptomatic   susceptible young_old
                         "antiviral",
                         time_start2, time_end2,
                         cvg_rate2, cvg_max_prop2);
    
    interv_vec.push_back(interv1);
    interv_vec.push_back(interv2);
    
    
    
    // ================================================================
    //     WORLD PARAMETERS
    // ================================================================
    
    // Area Units
    vector<ID>      id_au {0,1}; //2,3};
    vector<string>  name_au {"AUone", "AUtwo"}; //, "AUthree", "AUfour"};
    unsigned long   n_au = id_au.size();
    ID id_region = 0;
    string regionName = "RegionOne";
    vector<areaUnit> auvec = create_area_unit(id_au, name_au, id_region, regionName);
    
    uint mult = 1;
    
    
    // ====== PARAMETERS FOR STOCHASTIC WORLD =======
    
    
    // Vector of size distributions
    // (vector size is the number of AU,
    //  bc each AU has its own size distribution)
    vector<discrete_prob_dist<uint> > D_size_hh_vec;
    vector<discrete_prob_dist<uint> > D_size_wrk_vec;
    vector<discrete_prob_dist<uint> > D_size_school_vec;
    vector<discrete_prob_dist<uint> > D_size_pubt_vec;
    vector<discrete_prob_dist<uint> > D_size_hosp_vec;
    vector<discrete_prob_dist<uint> > D_size_other_vec;
    
    // Households sizes
    vector<uint> hh_size {1,2,3};
    vector<double> hh_size_proba {0.2, 0.5, 0.30};
    discrete_prob_dist<uint> D_size_hh(hh_size, hh_size_proba);
    
    // in this test file, assume both AU have same size distributions:
    D_size_hh_vec.push_back(D_size_hh);
    D_size_hh_vec.push_back(D_size_hh);
    
    // Workplace sizes
    vector<uint> wrk_size {5,20,40,60};
    vector<double> wrk_size_proba {0.55, 0.3, 0.1, 0.05};
    discrete_prob_dist<uint> D_size_wrk(wrk_size, wrk_size_proba);
    // in this test file, assume both AU have same size distributions:
    D_size_wrk_vec.push_back(D_size_wrk);
    D_size_wrk_vec.push_back(D_size_wrk);
    
    // Public transport sizes
    vector<uint> pubt_size {30,60,120};
    vector<double> pubt_size_proba {0.3,0.4,0.3};
    discrete_prob_dist<uint> D_size_pubt(pubt_size, pubt_size_proba);
    // in this test file, assume both AU have same size distributions:
    D_size_pubt_vec.push_back(D_size_pubt);
    D_size_pubt_vec.push_back(D_size_pubt);
    
    // School sizes
    vector<uint> school_size {100, 200, 300};
    vector<double> school_size_proba {0.7,0.2,0.1};
    discrete_prob_dist<uint> D_size_school(school_size, school_size_proba);
    // in this test file, assume both AU have same size distributions:
    D_size_school_vec.push_back(D_size_school);
    D_size_school_vec.push_back(D_size_school);
    
    // Hospital sizes
    vector<uint> hosp_size {50000};
    vector<double> hosp_size_proba {1};
    discrete_prob_dist<uint> D_size_hosp(hosp_size, hosp_size_proba);
    // in this test file, assume both AU have same size distributions:
    D_size_hosp_vec.push_back(D_size_hosp);
    D_size_hosp_vec.push_back(D_size_hosp);
    
    // other public placed sizes
    vector<uint> other_size {5,6,7};
    vector<double> other_size_proba {0.4,0.4,0.2};
    discrete_prob_dist<uint> D_size_other(other_size, other_size_proba);
    // in this test file, assume both AU have same size distributions:
    D_size_other_vec.push_back(D_size_other);
    D_size_other_vec.push_back(D_size_other);

    
    
    // number of each social place type
    vector<uint> n_hh       {1200*mult, 1001*mult};
    // WARNING:
    // make sure there are enough social places
    // of type different from 'household',
    // else some individual will not have destinations
    // and the code will crash.
    vector<uint> n_wrk      {300*mult, 301*mult};
    vector<uint> n_pubt     {400*mult, 301*mult};
    vector<uint> n_school   {200*mult, 201*mult};
    vector<uint> n_hosp     {1,1};
    vector<uint> n_other    {1000*mult, 501*mult};
    
    
    // ====== PARAMETERS FOR DETERMINSTIC WORLD =======
    
    
    // Vector of size distributions
    // (vector size is the number of AU,
    //  bc each AU has its own size distribution)
    vector<vector<uint> > size_hh_vec(n_au);
    vector<vector<uint> > size_wrk_vec(n_au);
    vector<vector<uint> > size_school_vec(n_au);
    vector<vector<uint> > size_pubt_vec(n_au);
    vector<vector<uint> > size_hosp_vec(n_au);
    vector<vector<uint> > size_other_vec(n_au);
    
    // * * NOTE * *
    // The goal is to provide a _stable_ definition of the world, once it has been identified.
    // (whereas the 'stochastic' construction is repeated - from scratch - every time,
    // generating potentially not consistent world)
    // Although the world construction is supposed to be determinsitic,
    // here I generate, with distributions because doing it manually would be tedious.
    
    vector<uint>   hh_val  = {1,2,3};
    vector<double> hh_pr   = {0.3,0.4,0.3};
   
    vector<uint>   wrk_val  = {5,25,50,100,};
    vector<double> wrk_pr   = {0.4,0.3,0.2,0.1};
    
    vector<uint>   school_val  = {100, 250, 500, 1000, 1500};
    vector<double> school_pr   = {0.2, 0.3, 0.3, 0.1, 0.1};
    
    vector<uint> pubt_val  = {10,50,100};
    vector<double> pubt_pr = {0.3,0.4,0.3};

    vector<uint> other_val  = {10,50,100, 500};
    vector<double> other_pr = {0.3,0.4,0.2, 0.1};

    discrete_prob_dist<uint> p_hh(hh_val,hh_pr);
    discrete_prob_dist<uint> p_wrk(wrk_val,wrk_pr);
    discrete_prob_dist<uint> p_school(school_val,school_pr);
    discrete_prob_dist<uint> p_pubt(pubt_val,pubt_pr);
    discrete_prob_dist<uint> p_other(other_val,other_pr);
    
    for(int a=0; a<n_au; a++) {
        size_hh_vec[a]     = p_hh.sample(1200*mult, a+1234);
        size_wrk_vec[a]    = p_wrk.sample(400*mult, a+1234);
        size_school_vec[a] = p_school.sample(1000*mult, a+1234);
        size_pubt_vec[a]   = p_pubt.sample(1200*mult, a+1234);
        size_other_vec[a]  = p_other.sample(1000*mult, a+1234);
        size_hosp_vec[a]   = {50000};
    }
    
    
    // Age distribution inside households (wether stoch or det world construction)
    
    vector< vector<discrete_prob_dist<uint> > > pr_age_hh;
    
    vector<uint> age_adult = {22, 33, 44, 55};
    vector<uint> age_child = {5, 11};
    vector<uint> age_all = age_child;
    age_all.insert(age_all.end(), age_adult.begin(),age_adult.end());
    
    vector<double> p_age_adult = {0.2, 0.5, 0.2, 0.1};
    vector<double> p_age_all_1 = {0.1, 0.1, 0.2, 0.3, 0.2, 0.1};
    vector<double> p_age_all_2 = {0.4, 0.3, 0.1, 0.1, 0.05, 0.05};
    
    discrete_prob_dist<uint> pr_age_hh_00(age_adult, p_age_adult);
    discrete_prob_dist<uint> pr_age_hh_10(age_adult, p_age_adult);
    discrete_prob_dist<uint> pr_age_hh_11(age_all,   p_age_all_1);
    discrete_prob_dist<uint> pr_age_hh_20(age_adult, p_age_adult);
    discrete_prob_dist<uint> pr_age_hh_21(age_adult, p_age_adult);
    discrete_prob_dist<uint> pr_age_hh_22(age_all,   p_age_all_2);
    
    vector<discrete_prob_dist<uint> > tmp;
    
    tmp.push_back(pr_age_hh_00);
    pr_age_hh.push_back(tmp);
    
    tmp.clear();
    tmp = {pr_age_hh_10,pr_age_hh_11};
    pr_age_hh.push_back(tmp);
    
    tmp.clear();
    tmp = {pr_age_hh_20,pr_age_hh_21,pr_age_hh_22};
    pr_age_hh.push_back(tmp);
    
    float unemployed_prop = 0.10;
    
    // === Schedules definition ===
    
    // Define the time slices
    // must be same for all schedules and must sum up to 1.0
    vector<double> timeslice {1.0/24, 4.0/24, 4.0/24, 1.0/24, 2.0/24, 12.0/24};
    
    vector<vector<string> > sched_des;
    vector<string>          sched_nam;
    
    sched_des.push_back({"SP_pubTransp", "SP_workplace", "SP_workplace", "SP_pubTransp", "SP_other", "SP_household"});
    sched_des.push_back({"SP_pubTransp", "SP_school",    "SP_school",    "SP_pubTransp", "SP_other", "SP_household"});
    sched_des.push_back({"SP_household", "SP_other",     "SP_other",     "SP_other",     "SP_other", "SP_household"});
    
    sched_nam.push_back("worker_sed");
    sched_nam.push_back("student");
    sched_nam.push_back("unemployed");
    
    vector<schedule> sched = build_all_schedules(sched_des, sched_nam, timeslice);
    
    
    // ================================================================
    //     RUN SIMULATION
    // ================================================================
    
    bool build_world_only           = false;
    bool stoch_world_construction   = false;
    Simulator sim ;
    
    if(stoch_world_construction){
        sim = run_stochWorld(auvec,
                             D_size_hh_vec,
                             D_size_wrk_vec,
                             D_size_pubt_vec,
                             D_size_school_vec,
                             D_size_hosp_vec,
                             D_size_other_vec,
                             pr_age_hh,
                             n_hh ,
                             n_wrk,
                             n_pubt ,
                             n_school,
                             n_hosp,
                             n_other,
                             unemployed_prop,
                             sched ,
                             MP,
                             start_time,
                             horizon,
                             i0,
                             interv_vec,
                             build_world_only);
    }
    else{
        sim = run_detWorld(auvec,
                           size_hh_vec,
                           size_wrk_vec,
                           size_pubt_vec,
                           size_school_vec,
                           size_hosp_vec,
                           size_other_vec,
                           pr_age_hh,
                           unemployed_prop,
                           sched ,
                           MP,
                           start_time,
                           horizon,
                           i0,
                           interv_vec,
                           build_world_only);
    }
    
    // ================================================================
    //     MANIPULATE RESULTS
    // ================================================================
    
    // insert code here ...
    double dummy = 1; dummy++;
    
    dcDataFrame x = export_world(sim.get_world());
    
    //x.display();
    //displayVector(census_ages(sim.get_world()));
}


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
                         bool build_world_only){
    
    
    stopif(auvec.size() != n_hh.size() ||
           auvec.size() != n_wrk.size() ||
           auvec.size() != n_pubt.size() ||
           auvec.size() != n_school.size() ||
           auvec.size() != n_hosp.size() ||
           auvec.size() != n_other.size(),
           "Inconsistent inputs.");
    
    Simulator sim;
    sim.set_modelParam(MP);
    
    sim.create_world(auvec,
                     D_size_hh, pr_age_hh,
                     D_size_wrk,
                     D_size_pubt,
                     D_size_school,
                     D_size_hosp,
                     D_size_other,
                     n_hh,
                     n_wrk,
                     n_pubt,
                     n_school,
                     n_hosp,
                     n_other,
                     unemployed_prop,
                     sched);
    
    sim.set_horizon(horizon);
    
    // Define the disease
    double dol_mean = MP.get_prm_double("dol_mean");
    double doi_mean = MP.get_prm_double("doi_mean");
    double doh_mean = MP.get_prm_double("doh_mean");
    double dol_var  = MP.get_prm_double("dol_var");
    double doi_var  = MP.get_prm_double("doi_var");
    double doh_var  = MP.get_prm_double("doh_var");
    
//    bool debug_mode = MP.get_prm_bool("debug_mode");
    
    disease flu("Influenza",
                dol_mean, doi_mean, doh_mean,
                dol_var, doi_var, doh_var);
    sim.set_disease(flu);
    
    // Interventions:
    for(int i=0; i<interv.size(); i++)
        sim.add_intervention(interv[i]);
    
    // Run the simulation:
    sim.set_start_time(start_time);
    sim.set_initial_prevalence(i0);
    
    sim.display_summary_info();
    
    if(!build_world_only) sim.run();
    else{
        cout << " World is built only, simulations were not requested."<<endl;
    }
    
    return sim;
}


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
                       bool build_world_only){
    
    Simulator sim;
    sim.set_modelParam(MP);
    
    sim.create_world_det(auvec,
                         size_hh, pr_age_hh,
                         size_wrk,
                         size_pubt,
                         size_school,
                         size_hosp,
                         size_other,
                         unemployed_prop,
                         sched);
    
    sim.set_horizon(horizon);
    
    // Define the disease
    double dol_mean = MP.get_prm_double("dol_mean");
    double doi_mean = MP.get_prm_double("doi_mean");
    double doh_mean = MP.get_prm_double("doh_mean");
    double dol_var  = MP.get_prm_double("dol_var");
    double doi_var  = MP.get_prm_double("doi_var");
    double doh_var  = MP.get_prm_double("doh_var");
    
    //    bool debug_mode = MP.get_prm_bool("debug_mode");
    
    disease flu("Influenza",
                dol_mean, doi_mean, doh_mean,
                dol_var, doi_var, doh_var);
    sim.set_disease(flu);
    
    // Interventions:
    for(int i=0; i<interv.size(); i++)
        sim.add_intervention(interv[i]);
    
    // Run the simulation:
    sim.set_start_time(start_time);
    sim.set_initial_prevalence(i0);
    
    sim.display_summary_info();
    
    if(!build_world_only) sim.run();
    else{
        cout << endl << " -- World is built only, simulations were not requested --"<<endl;
    }
    return sim;
}










Simulator test_naf(modelParam MP,
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
    
    Simulator sim;
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
    
    MP.add_prm_double ("vax_imm_hum_incr", 0.2);
    MP.add_prm_double ("vax_imm_cell_incr", 0.4);
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
    
    Simulator sim1 = test_naf(MP, horizon, n_indiv, i0,interv_test);
    
    dcDataFrame pop_final = sim1.get_world()[2].export_dcDataFrame();
    pop_final.display();
}







Simulator test_SEIR_vs_ODE(modelParam MP,
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
    
    Simulator sim;
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
    Simulator sim1 = test_SEIR_vs_ODE(MP,horizon,n_indiv,i0);
    
    sim1.get_world()[0].export_dcDataFrame().display();
    sim1.timeseries().display();
    
}



Simulator test_move_2_sp(modelParam MP,
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
    
    Simulator sim;
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
    Simulator sim1 = test_move_2_sp(MP, horizon, n_indiv, i0);
}




Simulator test_hospitalization(modelParam MP,
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
    
    Simulator sim;
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
    Simulator sim1 = test_hospitalization(MP, horizon, n_indiv, i0);
    
    dcDataFrame pop_final = sim1.get_world()[2].export_dcDataFrame();
    pop_final.display();
}



void test_move_transmission(){
    
    auto t0 = std::chrono::system_clock::now();
    
    // Build associated simulation:
    double horizon = 20;
    double sizereduction = 0.002; // Scale down world size compared to real world one
    
    Simulator sim;
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











