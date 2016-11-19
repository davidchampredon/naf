//
//  simulation.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-06.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "simulator.h"
#include "shapeFct.h"

void Simulator::base_constructor(){
    
    _start_time = -9999.99;
    
    _n_E  = 0;
    _n_Ia = 0;
    _n_Is = 0;
    _n_R  = 0;
    _n_H  = 0;
    _n_D  = 0;
    _incidence = 0;
    _n_treated.clear();
    _n_vaccinated.clear();
    _horizon = -999;
    
    _ts_times.clear();
    _ts_E.clear();
    _ts_Ia.clear();
    _ts_Is.clear();
    _ts_R.clear();
    _ts_D.clear();
    _ts_census_sp_time.clear();
    _ts_census_sp_id.clear();
    _ts_census_sp_type.clear();
    _ts_census_sp_nS.clear();
    _ts_census_sp_nE.clear();
    
    _ts_n_treated.clear();
    _ts_n_vaccinated.clear();
    
    _intervention.clear();
    
    _track_n_contacts.clear();
    _track_n_contacts_time.clear();
    _track_n_contacts_uid.clear();
    
    _wiw_ages.clear();
    _wiw_ages.resize(2);
}


Simulator::Simulator(){
    base_constructor();
}


void Simulator::create_world(vector<areaUnit> AU,
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
                             vector<uint> n_other,
                             float unemployed_prop,
                             vector<schedule> sched){
    
    cout << endl << " Creating world stochastically ... " << endl;
    
    // Create social places:
    
    world W = build_world(AU,
                          D_size_hh, pr_age_hh,
                          D_size_wrk,
                          D_size_pubt,
                          D_size_school,
                          D_size_hosp,
                          D_size_other,
                          n_hh, n_wrk, n_pubt, n_school, n_hosp, n_other);
    
    set_world(W);
    set_sp_other();
    
    // Define individuals' features:
    
    string dol_distrib = _modelParam.get_prm_string("dol_distrib");
    string doi_distrib = _modelParam.get_prm_string("doi_distrib");
    string doh_distrib = _modelParam.get_prm_string("doh_distrib");
    assign_dox_distribution(dol_distrib,
                            doi_distrib,
                            doh_distrib);
    
    assign_immunity_hum();
    assign_immunity_cell();
    assign_frailty();

    assign_schedules(_world, sched, unemployed_prop);
    
    check_schedules_consistency();
    
    cout << "... world created. " << endl;
}




void Simulator::create_world_det(vector<areaUnit> AU,
                                 vector<vector<uint> > size_hh,     // Households sizes
                                 vector< vector<discrete_prob_dist<uint> > > pr_age_hh,  // Age distribution inside households
                                 
                                 vector<vector<uint> > size_wrk,
                                 vector<vector<uint> > size_pubt,
                                 vector<vector<uint> > size_school,
                                 vector<vector<uint> > size_hosp,
                                 vector<vector<uint> > size_other,
                                 
                                 float unemployed_prop,
                                 vector<schedule> sched){
    
    cout << endl << " Creating world deterministically... " << endl;
    
    // Create social places:
    
    world W = build_world_det(AU,
                              size_hh, pr_age_hh,
                              size_wrk,
                              size_pubt,
                              size_school,
                              size_hosp,
                              size_other);
    
    set_world(W);
    set_sp_other();
    
    // Define individuals' features:
    
    string dol_distrib = _modelParam.get_prm_string("dol_distrib");
    string doi_distrib = _modelParam.get_prm_string("doi_distrib");
    string doh_distrib = _modelParam.get_prm_string("doh_distrib");
    assign_dox_distribution(dol_distrib,
                            doi_distrib,
                            doh_distrib);
    
    assign_immunity_hum();
    assign_immunity_cell();
    assign_frailty();
    
    assign_schedules(_world, sched, unemployed_prop);
    
    check_schedules_consistency();
    
    cout << "... world created. " << endl;
}



void Simulator::check_schedules_consistency(){
    /// Check that all individuals who have a
    /// given SP type in their schedule, are
    /// indeed linked to a SP of this type.
    
    vector<uint> cnt(SP_MAX, 0);
    
    for(uint k=0; k<_world.size(); k++){
        uint nk = (uint)_world[k].get_size();
        
        for(uint i=0; i<nk; i++){
            schedule sched =_world[k].get_indiv(i).get_schedule();
            uint n = (uint) sched.get_timeslice().size();
            for(uint s=0; s<n; s++){
                if(_world[k].get_indiv(i).find_dest(s,false) == __UNDEFINED_ID ){
                    SPtype q = sched.get_sp_type(s);
                    cnt[q]++;
                }
            }
        }
    }
    
    // display checks results:
    cout << endl <<  " === Schedules consistency === " <<endl<<endl;
    cout << " Number of individuals without\n a link to relevant SP type\n in its schedule (should be 0):  " <<endl<<endl;
    for(uint i=0; i<cnt.size(); i++){
        tabcout(SPtype2string(int2SPtype(i)), cnt[i],22);
    }
    cout<<endl;
    //stopif(sumElements(cnt) > 0, "Some individuals have no link to relevant SP.");
}




void Simulator::build_single_world(uint n_indiv){
    /// Build a very simple world:
    /// one social place with 'n_indiv' individuals.
    /// USED FOR TEST
    
    cout << endl << "Building single test world..."<<endl;
    
    string region1 = "Halton";
    ID id_region1 = 1;
    
    areaUnit A1(1, "singletown", id_region1, region1);
    
    vector<areaUnit> A {A1};
    
    // Schedule
    int nt = _modelParam.get_prm_uint("nt");
    vector<double> timeslice(nt, 1.0/nt); // must sum up to 1.0
    vector<SPtype> single_sed  {SP_household};
    schedule sched_worker_sed(single_sed, timeslice, "worker_sed");
    
    vector<schedule> sched {
        sched_worker_sed
    };
    
    // Disease stage durations:
    string dol_distrib = "exp";
    string doi_distrib = "exp";
    string doh_distrib = "exp";
    
    // individuals
    uint num_indiv			= (uint)(n_indiv * 3); // <-- make sure it's large enough
    vector<individual> many_indiv	= build_individuals(num_indiv,
                                                        sched,
                                                        dol_distrib,
                                                        doi_distrib,
                                                        doh_distrib );
    
    // type of social places
    vector<SPtype> spt {SP_household};
    
    // populate social places with individuals (built above)
    uint num_hh = 1;
    
    cout << "Number of household places: " << num_hh <<endl;
    
    // W A R N I N G
    // same order as type of social places definition ('spt')
    vector<uint> num_sp {num_hh};
    
    // Distribution of the size of each social place type:
    discrete_prob_dist<uint> p_hh({n_indiv},{1.0});
    
    vector< discrete_prob_dist<uint> > p_size {p_hh};
    
    // build the world:
    vector<socialPlace> W = build_world_simple(spt, num_sp, p_size, many_indiv, A);
    set_world(W);
    
    // initial population: Move everyone to its household!
    unsigned long N = _world.size();
    for (int k=0; k<N; k++)
    {
        if(_world[k].get_type()==SP_household)
        {
            uint n_linked_k = _world[k].n_linked_indiv();
            vector<ID> linked_ids   = _world[k].get_linked_indiv_id();
            
            stopif(linked_ids.size() != n_linked_k, "Book keeping problem with linked IDs");
            
            for (uint i=0; i<n_linked_k; i++)
            {
                //ID curr_indiv_ID = _world[k].get_linked_indiv_id()[i];
                ID curr_indiv_ID = linked_ids[i];
                individual tmp = get_indiv_with_ID(curr_indiv_ID, many_indiv);
                
                bool debugcode = false;
                if(debugcode){
                    // Checks (remove for better speed)
                    ID id_hh = tmp.get_id_sp_household();
                    stopif(id_hh == __UNDEFINED_ID, "at least one individual has no linked household!");
                    stopif(id_hh != k, "Not consistent linkage!");
                    _world[id_hh].add_indiv(tmp);
                    // -----
                }
                
                // Faster version:
                if(!debugcode) _world[k].add_indiv(tmp);
            }
        }
    }
    cout << "... test world built." << endl;
}


void Simulator::build_test_2_sp(uint n_indiv){
    /// Build a world with only 2 social places
    /// to test movements
    
    cout << endl << "Building test world 2 SP..."<<endl;
    
    bool debugcode = _modelParam.get_prm_bool("debug_mode");
    
    string region1 = "Halton";
    ID id_region1 = 1;
    
    areaUnit A1(1, "Oakville", id_region1, region1);
    areaUnit A2(2, "Burlington", id_region1, region1);
    
    vector<areaUnit> A {A1,A2};
    
    // Schedule
    vector<double> timeslice {12.0/24, 12.0/24}; // must sum up to 1.0
    vector<SPtype> worker_sed  {SP_workplace, SP_household};
    
    schedule sched_worker_sed(worker_sed, timeslice, "worker_sed");
    
    vector<schedule> sched {
        sched_worker_sed,
    };
    
    // Disease stage durations:
    string dol_distrib = "exp";
    string doi_distrib = "exp";
    string doh_distrib = "exp";
    // individuals
    uint num_indiv			= (uint)(n_indiv*2); // <-- make sure it's large enough
    vector<individual> many_indiv	= build_individuals(num_indiv,
                                                        sched,
                                                        dol_distrib,
                                                        doi_distrib,
                                                        doh_distrib);
    
    // type of social places
    vector<SPtype> spt {
        SP_workplace,
        SP_household
    };
    
    // populate social places with individuals (built above)
    uint num_hh      = (uint)(1 );
    uint num_biz     = (uint)(1 );
    
    cout << "Number of business places:  " << num_biz <<endl;
    cout << "Number of household places: " << num_hh <<endl;
    
    // Number of social places, by type.
    // * *   W A R N I N G   * *
    // same order as type of social places definition ('spt')
    vector<uint> num_sp {
        num_biz,
        num_hh
    };
    
    // Distribution of the size of each social place type:
    discrete_prob_dist<uint> p_workPlace({n_indiv/2},  {1.0});
    discrete_prob_dist<uint> p_hh({n_indiv},  {1.0});
    
    vector<discrete_prob_dist<uint> > p_size {
        p_workPlace,
        p_hh
    };
    
    // build the world:
    vector<socialPlace> W = build_world_simple(spt,
                                               num_sp,
                                               p_size,
                                               many_indiv, A);
    set_world(W);
    
    // initial population: Move everyone to its household!
    
    unsigned long N = _world.size();
    for (int k=0; k<N; k++)
    {
        if(_world[k].get_type()==SP_household)
        {
            uint n_linked_k = _world[k].n_linked_indiv();
            for (uint i=0; i<n_linked_k; i++)
            {
                ID curr_indiv_ID = _world[k].get_linked_indiv_id()[i];
                individual tmp = get_indiv_with_ID(curr_indiv_ID, many_indiv);
                
                
                if(debugcode){
                    // Checks (remove for better speed)
                    ID id_hh = tmp.get_id_sp_household();
                    stopif(id_hh == __UNDEFINED_ID, "at least one individual has no linked household!");
                    stopif(id_hh != k, "Not consistent linkage!");
                    _world[id_hh].add_indiv(tmp);
                    // -----
                }
                
                // Faster version:
                if(!debugcode) _world[k].add_indiv(tmp);
            }
        }
    }
    cout << "... test world built."<<endl;
    
}



void Simulator::build_test_hospitalization(uint n_indiv){
    /// Build a world
    /// to test hospitalization and movements
    
    cout << endl << "Building test world hospitalization..."<<endl;
    
    bool debugcode = _modelParam.get_prm_bool("debug_mode");
    
    string region1 = "Halton";
    ID id_region1 = 1;
    
    areaUnit A1(1, "Oakville", id_region1, region1);
    areaUnit A2(2, "Burlington", id_region1, region1);
    areaUnit A3(2, "Hospital", id_region1, region1);
    
    vector<areaUnit> A {A1,A2,A3};
    
    // type of social places
    // existing in the simulated world:
    vector<SPtype> spt {
        SP_workplace,
        SP_household,
        SP_hospital
    };
    
    
    // === Schedules ===
    
    // Define the time slices
    // must be same for all schedules and must sum up to 1.0
    vector<double> timeslice {6.0/24, 18.0/24};
    
    // type of schedules:
    vector<SPtype> worker_sed  {SP_workplace, SP_household};
    vector<SPtype> stayhome  {SP_household, SP_household};
    
    schedule sched_worker_sed(worker_sed, timeslice, "worker_sed");
    schedule sched_stayhome(stayhome, timeslice, "stayhome");
    
    // Schedules used in the simulation:
    vector<schedule> sched {
        sched_worker_sed,
        sched_stayhome
    };
    
    
    // Disease stage durations:
    string dol_distrib = _modelParam.get_prm_string("dol_distrib");
    string doi_distrib = _modelParam.get_prm_string("doi_distrib");
    string doh_distrib = _modelParam.get_prm_string("doh_distrib");
    
    // individuals
    uint num_indiv = (uint)(n_indiv*2); // <-- make sure it's large enough
    vector<individual> many_indiv = build_individuals(num_indiv,
                                                      sched,
                                                      dol_distrib,
                                                      doi_distrib,
                                                      doh_distrib);
    
    
    // populate social places with individuals (built above)
    uint num_hh      = (uint)(1 );
    uint num_biz     = (uint)(1 );
    uint num_hosp    = (uint)(1 );
    
    cout << "Number of business places:  " << num_biz << endl;
    cout << "Number of household places: " << num_hh << endl;
    cout << "Number of hospital places: " << num_hosp << endl;
    
    // Number of social places, by type.
    // * *   W A R N I N G   * *
    // same order as type of social places definition ('spt')
    vector<uint> num_sp {
        num_biz,
        num_hh,
        num_hosp
    };
    
    // Distribution of the size of each social place type:
    discrete_prob_dist<uint> p_workPlace({n_indiv/2},  {1.0});
    discrete_prob_dist<uint> p_hh({n_indiv},  {1.0});
    discrete_prob_dist<uint> p_hosp({0},  {1.0});
    
    vector<discrete_prob_dist<uint> > p_size {
        p_workPlace,
        p_hh,
        p_hosp
    };
    
    // build the world:
    vector<socialPlace> W = build_world_simple_2(
                                                 many_indiv,
                                                 A,
                                                 sched);
    set_world(W);
    
    // initial population: Move everyone to its household!
    
    unsigned long N = _world.size();
    for (int k=0; k<N; k++)
    {
        if(_world[k].get_type()==SP_household)
        {
            uint n_linked_k = _world[k].n_linked_indiv();
            for (uint i=0; i<n_linked_k; i++)
            {
                ID curr_indiv_ID = _world[k].get_linked_indiv_id(i);
                individual tmp = get_indiv_with_ID(curr_indiv_ID, many_indiv);
                
                
                if(debugcode){
                    // Checks (remove for better speed)
                    ID id_hh = tmp.get_id_sp_household();
                    stopif(id_hh == __UNDEFINED_ID, "at least one individual has no linked household!");
                    stopif(id_hh != k, "Not consistent linkage!");
                    _world[id_hh].add_indiv(tmp);
                    // -----
                }
                
                // Faster version:
                if(!debugcode) _world[k].add_indiv(tmp);
            }
        }
    }
    cout << "... test world built."<<endl;
    
}




void Simulator::build_test(uint n_indiv){
    /// Build a world for test
    
    cout << endl << "Building test world..."<<endl;
    
    bool debugcode = _modelParam.get_prm_bool("debug_mode");
    
    string region1 = "Halton";
    ID id_region1 = 1;
    
    areaUnit A1(1, "Oakville", id_region1, region1);
    areaUnit A2(2, "Burlington", id_region1, region1);
    areaUnit A3(2, "Hospital", id_region1, region1);
    
    vector<areaUnit> A {A1,A2,A3};
    
    // type of social places
    // existing in the simulated world:
    vector<SPtype> spt {
        SP_workplace,
        SP_household,
        SP_hospital
    };
    
    
    // === Schedules ===
    
    // Define the time slices
    // must be same for all schedules and must sum up to 1.0
    vector<double> timeslice {6.0/24, 18.0/24};
    
    // type of schedules:
    vector<SPtype> worker_sed  {SP_workplace, SP_household};
    vector<SPtype> stayhome  {SP_household, SP_household};
    
    schedule sched_worker_sed(worker_sed, timeslice, "worker_sed");
    schedule sched_stayhome(stayhome, timeslice, "stayhome");
    
    // Schedules used in the simulation:
    vector<schedule> sched {
        sched_worker_sed,
        sched_stayhome
    };
    
    
    // Disease stage durations:
    string dol_distrib = _modelParam.get_prm_string("dol_distrib");
    string doi_distrib = _modelParam.get_prm_string("doi_distrib");
    string doh_distrib = _modelParam.get_prm_string("doh_distrib");
    
    // individuals
    uint num_indiv = (uint)(n_indiv*2); // <-- make sure it's large enough
    vector<individual> many_indiv = build_individuals(num_indiv,
                                                      sched,
                                                      dol_distrib,
                                                      doi_distrib,
                                                      doh_distrib);
    
    
    // populate social places with individuals (built above)
    uint num_hh      = (uint)(1 );
    uint num_biz     = (uint)(1 );
    uint num_hosp    = (uint)(1 );
    
    cout << "Number of business places:  " << num_biz <<endl;
    cout << "Number of household places: " << num_hh <<endl;
    cout << "Number of hospital places: " << num_hosp <<endl;
    
    // Number of social places, by type.
    // * *   W A R N I N G   * *
    // same order as type of social places definition ('spt')
    vector<uint> num_sp {
        num_biz,
        num_hh,
        num_hosp
    };
    
    // Distribution of the size of each social place type:
    discrete_prob_dist<uint> p_workPlace({n_indiv/2},  {1.0});
    discrete_prob_dist<uint> p_hh({n_indiv},  {1.0});
    discrete_prob_dist<uint> p_hosp({0},  {1.0});
    
    vector<discrete_prob_dist<uint> > p_size {
        p_workPlace,
        p_hh,
        p_hosp
    };
    
    // build the world:
    vector<socialPlace> W = build_world_simple_2(
                                                 many_indiv,
                                                 A,
                                                 sched);
    set_world(W);
    
    // initial population: Move everyone to its household!
    
    unsigned long N = _world.size();
    for (int k=0; k<N; k++)
    {
        if(_world[k].get_type()==SP_household)
        {
            uint n_linked_k = _world[k].n_linked_indiv();
            for (uint i=0; i<n_linked_k; i++)
            {
                ID curr_indiv_ID = _world[k].get_linked_indiv_id(i);
                individual tmp = get_indiv_with_ID(curr_indiv_ID, many_indiv);
                
                
                if(debugcode){
                    // Checks (remove for better speed)
                    ID id_hh = tmp.get_id_sp_household();
                    stopif(id_hh == __UNDEFINED_ID, "at least one individual has no linked household!");
                    stopif(id_hh != k, "Not consistent linkage!");
                    _world[id_hh].add_indiv(tmp);
                    // -----
                }
                
                // Faster version:
                if(!debugcode) _world[k].add_indiv(tmp);
            }
        }
    }
    cout << "... test world built."<<endl;
    
}




void Simulator::build_test_world(double sizereduction){
    
    /// Build a test world with individuals LINKED and PRESENT to social places.
    /// NOTE: this function is for test, it will eventually be useless...
    
    cout << endl << "Building test world..."<<endl;
    
    string region1 = "Halton";
    ID id_region1 = 1;
    
    areaUnit A1(1, "Oakville", id_region1, region1);
    areaUnit A2(2, "Burlington", id_region1, region1);
    areaUnit A3(3, "Hamilton", id_region1, region1);
    
    vector<areaUnit> A {A1,A2,A3};
    
    // Schedule
    vector<double> timeslice {2.0/24, 8.0/24, 2.0/24, 12.0/24}; // must sum up to 1.0
    vector<SPtype> worker_sed  {SP_pubTransp, SP_workplace, SP_pubTransp, SP_household};
    vector<SPtype> student     {SP_pubTransp, SP_school,    SP_pubTransp, SP_household};
    
    schedule sched_worker_sed(worker_sed, timeslice, "worker_sed");
    schedule sched_student(student, timeslice, "student");
    
    //	vector<SPtype> worker_trav {SP_pubTransp, SP_workplace, SP_other, SP_household};
    //	schedule sched_worker_trav(worker_trav, timeslice, "worker_trav");
    
    vector<schedule> sched {
        sched_worker_sed,
        sched_student
    };
    
    // Disease stage durations:
    string dol_distrib = "exp";
    string doi_distrib = "exp";
    string doh_distrib = "exp";
    
    // individuals
    uint num_indiv			= (uint)(1e7 * sizereduction); // <-- make sure it's large enough
    vector<individual> many_indiv	= build_individuals(num_indiv,
                                                        sched,
                                                        dol_distrib,
                                                        doi_distrib,
                                                        doh_distrib);
    
    // type of social places
    vector<SPtype> spt {
        SP_pubTransp,
        SP_workplace,
        SP_household,
        SP_school
    };
    
    // populate social places with individuals (built above)
    uint num_pubTr   = (uint)(2000 * sizereduction);
    uint num_biz     = (uint)(180e3 * sizereduction);
    uint num_hh      = (uint)(1500e3 * sizereduction);
    uint num_school  = (uint)(2600 * sizereduction);
    
    cout << "Number of public transportations places: " << num_pubTr <<endl;
    cout << "Number of business places: " << num_biz <<endl;
    cout << "Number of household places: " << num_hh <<endl;
    cout << "Number of school places: " << num_school <<endl;
    
    // W A R N I N G
    // same order as type of social places definition ('spt')
    vector<uint> num_sp {
        num_pubTr,
        num_biz,
        num_hh,
        num_school
    };
    
    // Distribution of the size of each social place type:
    
    discrete_prob_dist<uint> p_pubTransp({20,30,60},{0.6,0.3,0.1});
    discrete_prob_dist<uint> p_workPlace({3,7,15,30,75,200},{0.6,0.15,0.15,0.07,0.02,0.01});
    discrete_prob_dist<uint> p_hh({1,2,3,4,5,6,7,8},{0.23, 0.34, 0.16, 0.15, 0.06, 0.03, 0.02, 0.01});
    discrete_prob_dist<uint> p_school({250,500,750,1000,1250,1500},{0.60,0.25,0.10,0.03, 0.01,0.01});
    
    vector<discrete_prob_dist<uint> > p_size {
        p_pubTransp,
        p_workPlace,
        p_hh,
        p_school
    };
    
    // build the world:
    vector<socialPlace> W = build_world_simple(spt, num_sp, p_size, many_indiv, A);
    set_world(W);
    
    // initial population: Move everyone to its household!
    unsigned long N = _world.size();
    for (int k=0; k<N; k++)
    {
        if(_world[k].get_type()==SP_household)
        {
            uint n_linked_k = _world[k].n_linked_indiv();
            for (uint i=0; i<n_linked_k; i++)
            {
                ID curr_indiv_ID = _world[k].get_linked_indiv_id()[i];
                individual tmp = get_indiv_with_ID(curr_indiv_ID, many_indiv);
                
                bool debugcode = false;
                if(debugcode){
                    // Checks (remove for better speed)
                    ID id_hh = tmp.get_id_sp_household();
                    stopif(id_hh == __UNDEFINED_ID, "at least one individual has no linked household!");
                    stopif(id_hh != k, "Not consistent linkage!");
                    _world[id_hh].add_indiv(tmp);
                    // -----
                }
                
                // Faster version:
                if(!debugcode) _world[k].add_indiv(tmp);
            }
        }
    }
    cout << "... test world built."<<endl;
}


void Simulator::time_update(double dt){
    /// Make all relevant updates when time is advanced by 'dt'
    
    // Main timer:
    _current_time += dt;
    
    // Update individuals' clock:
    for (uint k=0; k<_world.size(); k++) {
        _world[k].time_update(dt);
    }
    update_pop_count();
}


void Simulator::test(){
    
    //    _world[0]._indiv_S[0]->set_immunity(0.12345);
    
}

void Simulator::initial_infections(uint i0){
    
    vector<ID> id_pres = at_least_one_indiv_present(_world);

    vector<ID> sp_initially_infected;
    vector<uint> I0;
    for(uint m=0; m<i0; m++){
        // Randomly select where infections are seeded
        std::uniform_int_distribution<uint> unif(0,(uint)id_pres.size());
        uint idx_rnd = unif(_RANDOM_GENERATOR);
        sp_initially_infected.push_back(id_pres[idx_rnd]);
        I0.push_back(1);
    }
    define_all_id_tables();
    seed_infection(sp_initially_infected, I0);
}


void Simulator::run(){
    /// Run the simulated epidemic
    
    bool debug_mode = _modelParam.get_prm_bool("debug_mode");
    if(debug_mode) cout << endl << endl << " ======= START SIMULATION ======" <<endl<<endl;
    
    // TO DO: CHANGE THAT, IT's UGLY AND DANGEROUS
    ID ii = at_least_one_indiv_present(_world)[0];
    vector<double> timeslice = _world[ii].get_indiv()[0].get_schedule().get_timeslice();
    // - - - - - - - - -
    
    unsigned long nts = timeslice.size();
    uint k = 0;
    
    if(debug_mode) check_book_keeping();
    define_all_id_tables();
    
    // Interventions
    
    float treat_doi_reduc   = _modelParam.get_prm_double("treat_doi_reduc");
    float vax_imm_hum_incr  = _modelParam.get_prm_double("vax_imm_hum_incr");
    float vax_imm_cell_incr = _modelParam.get_prm_double("vax_imm_cell_incr");
    float vax_frail_incr    = _modelParam.get_prm_double("vax_frail_incr");
    float vax_lag           = _modelParam.get_prm_double("vax_lag_full_efficacy");
    
    // Display simulator's information
    // before running simulation:
    display_summary_info();
    
    // set-up before time loop:
    double p_move    = _modelParam.get_prm_double("proba_move");
    double red_sympt = _modelParam.get_prm_double("proba_move_reduc_sympt");
    define_contactAssort();
    count_targeted_by_intervention();
    
    // ----- MAIN LOOP FOR TIME ------
    
    _current_time = _start_time;
    double tiny_time = 1E-6;
    
    while (_current_time < _horizon )
    {
        // If zero prevalence after the epidemic has starte, then stop.
        if(!at_least_one_infected() && _current_time>tiny_time) {break;}
        
        // Epidemic starts at time t = 0.0:
        if (fabs(_current_time-0.0)<=tiny_time){
            initial_infections(_initial_prevalence);
        }
        
        uint idx_timeslice = k % nts;
        double dt = timeslice[idx_timeslice];
        
        if(debug_mode){
            cout.precision(3);
            cout << "iter = " << k;
            cout << " ; sched idx = " << idx_timeslice;
            cout << " (length=" << dt << ")";
            cout << " ; time = " << _current_time;
            cout << " ; prev = " << _prevalence;
            cout <<endl;
            cout.precision(15);
        }
        else{
            cout.precision(3);
            if(k % 1 == 0){
                cout << "simulation time: " << _current_time;
                cout << "  (prevalence = "<<_prevalence<<")"<< endl;
            }
            cout.precision(15);
        }
        
        update_ts_census_SP();
        if(debug_mode) cout << "update_ts_census_SP done." <<endl;
        
        discharge_hospital(idx_timeslice);
        define_all_id_tables();
        death_hospital();
        if(debug_mode) cout << "hospitalization stuff done." <<endl;
        
        if(p_move>0) {
            move_individuals_sched(idx_timeslice, p_move, red_sympt);
            if(debug_mode) cout << "movements done."<<endl;
            define_all_id_tables();
            if(debug_mode) {
                cout << "id tables updated done."<<endl;
                check_book_keeping();
                cout << "book keeping checks done."<<endl;
            }
        }
        
        transmission_world(dt);
        update_pop_count();
        if(debug_mode) cout << "transmissions done." <<endl;
        
        hospitalize();
        define_all_id_tables(); // <-- check if this is necessary here
        update_pop_count(); // <-- check if this is necessary here
        if(debug_mode) cout << "hospitalization(2) done." <<endl;
        
        change_rnd_sp_other();
        define_all_id_tables();
        if(debug_mode) cout << "change other SP done." <<endl;
       
        // interventions:
        
        for (uint k=0; k<_world.size(); k++) {
            activate_interventions(k, dt,
                                   treat_doi_reduc,
                                   vax_imm_hum_incr,
                                   vax_imm_cell_incr,
                                   vax_frail_incr,
                                   vax_lag);
        }
        
        if(debug_mode) cout << "interventions done." <<endl;
        
        if(at_least_one_vaccination_intervention())
            update_immunities();
        if(debug_mode) cout << "immunity update done." <<endl;
        
        // Record for time series:
        timeseries_update();
        
        // Advance time:
        time_update(dt);
        k++;
        
        if(debug_mode) {
            check_book_keeping();
            cout << "all actions for this time step done." << endl<<endl;
        }
    }// end-while
    
    if(debug_mode){
        cout << endl << endl << "Simulator completed."<< endl;
        cout << "time series of incidence:"<<endl;
        displayVector(_ts_incidence);
        
        cout <<"Final size: " << sumElements(_ts_incidence)<<endl;
        cout << endl << "time series of hospitalizations:"<<endl;
        displayVector(_ts_H);
        cout << endl << "time series of deaths:"<<endl;
        displayVector(_ts_D);
    }
    cout << endl;
}


void Simulator::set_world(world w){
    _world = w;
    update_pop_count();
}


void Simulator::set_disease(const disease &d){
    /// Set the disease 'd' to all individuals in all social places
    for (uint k=0; k<_world.size(); k++) {
        _world[k].set_disease_to_all_indiv(d);
    }
}


void Simulator::move_one_individual(uint pos_indiv, ID from, ID to){
    
    individual tmp = _world[from].get_indiv(pos_indiv);
    
    // add individual at destination
    _world[to].add_indiv(tmp);
    // remove this individual (in pos_indiv^th position in '_indiv' vector) from here
    _world[from].remove_indiv(pos_indiv);
    // DEBUG
    // ID tmpid =_world[from].get_indiv(pos_indiv).get_id();
    //cout << "DEBUG --> id_" <<tmpid << "moved from SP_"<<from << " to SP_"<<to <<endl;
}


void Simulator::move_individuals_sched(uint idx_timeslice,
                                       double proba,
                                       double red_sympt){
    bool debug_mode = _modelParam.get_prm_bool("debug_mode");
    unsigned long N = _world.size();
    std::uniform_real_distribution<double> unif(0.0,1.0);
    
    for (int k=0; k<N; k++)
    {
        //if(debug_mode) cout<<"moving indiv in SP "<<k<<"/"<< N <<" ... ";
        uint n = (uint)_world[k].get_size();
        //if(debug_mode && n==0) cout << " (empty)."<<endl;
        if(n>0){
            
            // * * * IMPORTANT * * *
            // Must run this loop
            // in _descending_ order, else
            // it messes up the pointers vector deletion
            //
            for (uint i=n; (i--) > 0; )
            {
                // Check if individual is hospitalized
                // (hospitalized cannot move to other places!)
                if(! _world[k].get_indiv(i).is_hosp())
                {
                    // Retrieve its actual destination
                    ID id_dest = _world[k].find_dest(i, idx_timeslice);

                    if( id_dest == __UNDEFINED_ID){
                        string errmsg = "Undefined destination for indiv ID_";
                        errmsg += to_string(_world[k].get_indiv(i).get_id());
                        errmsg += " who is currently in SP_id_";
                        errmsg += to_string(k);
                        _world[k].get_indiv(i).displayInfo();
                        stopif(id_dest==__UNDEFINED_ID,errmsg);
                    }
                    
                    if ( id_dest != k ){
                        // Draw the chance move will actually happen:
                        double u = unif(_RANDOM_GENERATOR);
                        
                        // Probability reduction if symptomatic:
                        double m = 1.0;
                        if(_world[k].get_indiv(i).is_symptomatic()) m = red_sympt;
                        
                        // Draw if move happens:
                        if ( u < proba * m ) move_one_individual(i, k, id_dest);
                    }
                } // end-if-not-hospitalized
            }
        }
        if(debug_mode) cout<<" done." << endl;
    } // end-for-k-socialPlace
} // end-function


void Simulator::move_individuals(const SPtype sptype, double proba){
    
    /// Move individuals across social places
    
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    
    for (int k=0; k<_world.size(); k++)
    {
        vector<uint> pos2move; // idx position to move (must be done at the end of the loop, bc vector keeps on changing size!)
        for (int i=0; i<_world[k].get_size(); i++)
        {
            // Draw the chance move will actually happen:
            double u = unif(_RANDOM_GENERATOR);
            if ( u < proba )
            {
                individual tmp = _world[k].get_indiv(i);
                
                ID id_dest = __UNDEFINED_ID;
                
                if(sptype == SP_household)	id_dest = tmp.get_id_sp_household();
                if(sptype == SP_workplace)	id_dest = tmp.get_id_sp_workplace();
                if(sptype == SP_school)	    id_dest = tmp.get_id_sp_school();
                if(sptype == SP_other)		id_dest = tmp.get_id_sp_other();
                if(sptype == SP_hospital)	id_dest = tmp.get_id_sp_hospital();
                if(sptype == SP_pubTransp)	id_dest = tmp.get_id_sp_pubTransp();
                
                if(id_dest != __UNDEFINED_ID &&
                   id_dest != k)
                {
                    // add indiv to destination
                    _world[id_dest].add_indiv(tmp);
                    // record its position for future deletion
                    pos2move.push_back(i);
                }
            }
        }
        // remove from this SP the individuals that moved:
        if(pos2move.size()>0) _world[k].remove_indiv(pos2move);
    }
}


vector<uint> Simulator::draw_n_contacts(uint k,
                                        double dt,
                                        string infectious_type){
    /// Draw the (random) number of contacts
    /// for all infectious individuals of a given type.
    
    uint n=0;
    if(infectious_type == "Is") n = (uint)_world[k]._indiv_Is.size();
    if(infectious_type == "Ia") n = (uint)_world[k]._indiv_Ia.size();
    
    vector<uint> n_contacts(n,0); // number of contacts for each infectious individual
    
    uint nS = _world[k].get_n_S();
    uint total_contacts = 0;
    
    for (uint i=0; (i<n) && (total_contacts <= nS) ; i++) {
        // draw the contact rate
        // based on the individual and the current social place:
        individual* indiv = nullptr;
        if(infectious_type == "Is") indiv = _world[k]._indiv_Is[i];
        if(infectious_type == "Ia") indiv = _world[k]._indiv_Ia[i];
        
        double cr = draw_contact_rate(indiv, k);
        
        // Draw the actual number of contacts
        std::poisson_distribution<> poiss(cr * dt);
        uint rnd_contact = poiss(_RANDOM_GENERATOR);
        
        // tracking (used for debug)
        _track_n_contacts.push_back(rnd_contact);
        _track_n_contacts_time.push_back(_current_time);
        if(infectious_type == "Is")       _track_n_contacts_uid.push_back(_world[k]._indiv_Is[i]->get_id());
        else if (infectious_type == "Ia") _track_n_contacts_uid.push_back(_world[k]._indiv_Ia[i]->get_id());
        
        //DEBUG:
//        cout<<"base_cr = "<<_modelParam.get_prm_double("contact_rate_mean");
//        cout << " ; cr = "<< cr << " ; cr*dt = "<<cr*dt << " ; nc = " << rnd_contact <<endl;
        
        
        // If too many contact drawn,
        // reduce the latest one drawn,
        // and exit this loop.
        // (hence, the following infectious will have ZERO contact)
        total_contacts += rnd_contact;
        if ( total_contacts > nS ){
            rnd_contact -= (total_contacts - nS);
        }
        n_contacts[i] = rnd_contact;
    }
    return n_contacts;
}



void Simulator::set_sp_other(){
        
    for (uint k=0; k<_world.size(); k++) {
        if(_world[k].get_type() == SP_other){
            _sp_other.push_back(&_world[k]);
        }
    }
}


void Simulator::set_sp_other_link(uint k, uint pos_indiv, socialPlace &sp){
    /// link social place 'other' whose ID is 'id_sp'
    /// with the individual in _position_ 'pos' currently in kth social place
    
    _world[k].set_sp_other_link(pos_indiv, sp);
}


socialPlace* Simulator::pick_rnd_sp_other(){
    /// Choose randomly a social place
    /// of type 'other'
    
    std::uniform_int_distribution<uint> unif(0,(uint)_sp_other.size()-1);
    uint rnd_idx = unif(_RANDOM_GENERATOR);
    
    return _sp_other[rnd_idx];
}


double Simulator::calc_proba_transmission(individual *infectious,
                                          individual *susceptible){
    
    // TO DO: implement something more elaborate!?
    
    // Susceptible side (acquisition risk):
    double p_susc = 1.0 - susceptible->get_immunity_hum();
    
    // infectious side (transmission risk):
    double p_inf = 1.0;
    if(! infectious->is_symptomatic())
        p_inf = _modelParam.get_prm_double("asymptom_infectiousness_ratio");
    
    // Infectiousness is reduced
    // if patient is treated:
    if (infectious->is_treated()){
        double m = _modelParam.get_prm_double("treat_reduc_infect_mean");
        if(m>0){
            double reduc = beta_distribution(1.0, 1/m - 1.0, _RANDOM_GENERATOR);
            p_inf = p_inf * (1-reduc);
        }
    }
    return p_susc * p_inf;
}



double Simulator::calc_proba_symptomatic(float immunity, float frailty){
    /// Probability to be symptomatic given
    /// an individual's immunity and frailty
    
    // TO DO: more sophisticated!
    double res = frailty * (1-immunity); //(frailty * (1-immunity) < 0.5)? 0.2 : 1.0;
    return res;
}


double Simulator::calc_proba_hospitalized(float frailty)
{
    return frailty * _modelParam.get_prm_double("proba_hosp") ;
}


double Simulator::calc_proba_death(float frailty){
    // TO DO: more sophisticated!
    double p        = _modelParam.get_prm_double("proba_death_prm_1");
    double thres    = _modelParam.get_prm_double("proba_death_prm_2");
    double p_hi     = _modelParam.get_prm_double("proba_death_prm_3");
    if (frailty > thres) p = p_hi;
    return p;
}



vector< vector<uint> > Simulator::draw_contacted_S(uint k,
                                                   vector<uint> n_contacts,
                                                   string infectious_type){
    /// Randomly select the positions in '_indiv_S'
    /// for susceptibles that will be
    /// contacted by infectious individuals of a
    /// given infectious type (e.g. symptomatic or not).
    
    /// The returned format is a vector of vector, where
    /// column i can be seen as the list of susceptible
    /// contacted by the ith infector:
    ///
    /// infector_1     infector_2     infector_3    ...
    /// --------------------------------------------------
    /// indiv_S[.]     indiv_S[.]     indiv_S[.]    ...
    /// indiv_S[.]                    indiv_S[.]    ...
    /// indiv_S[.]                    indiv_S[.]    ...
    /// indiv_S[.]                                  ...
    
    uint n_inf = 0;
    if(infectious_type == "Is")         n_inf = (uint)_world[k]._indiv_Is.size();
    else if(infectious_type == "Ia")    n_inf = (uint)_world[k]._indiv_Ia.size();
    
    uint n_susc = (uint)_world[k]._indiv_S.size();
    
    stopif(sumElements(n_contacts)>n_susc, "Asking to draw more susceptibles than it exists!");
    
    vector< vector<uint> > selected_S(n_inf);
    
    // For all infectious indiv,
    // sample (with replacement)
    // which susceptibles will be contacted:
    std::uniform_int_distribution<uint> unif(0,n_susc-1);
    
    for(uint i=0; i<n_inf; i++){
        for(uint j=0; j<n_contacts[i]; j++)
        {
            uint rnd_pos = unif(_RANDOM_GENERATOR);
            selected_S[i].push_back(rnd_pos);
        }
    }
    return selected_S;
}


double age_contact_elem(double x,double y,
                        double a,double b,
                        double w,double q,double r){
    
    double tmp1 = pow( (x-a)-(y-b)-q ,2)/w/w;
    double tmp2 = (pow(x-a,2) + pow(y-b,2))/r/r ;
    double res  = exp(-tmp1 -tmp2);
    return res;
}


void Simulator::define_contactAssort(){
    /// Defines the age assortativity for contacts
    
    _contactAssort.resize(101);
    
    for(uint i=0; i<=100; i++){
        for(uint j=0; j<=i; j++){
            double tmp = 0.0;
            
            tmp += age_contact_elem(i,j, 0,  0,  9,  28, 45);    // parents/children
            tmp += age_contact_elem(i,j, 10, 10, 9,  0,  28);    // children/children
            tmp += age_contact_elem(i,j, 45, 45, 20, 0,  17);    // adults/adults
            tmp += age_contact_elem(i,j, 75, 75, 15, 0,  15);    // seniors
            tmp += age_contact_elem(i,j, 0,  0,  10, 60, 55);    // seniors/children

            // Symetrical matrix:
            _contactAssort(i,j) = tmp;
            _contactAssort(j,i) = tmp;
        }
    }
    // Normalize:
    double z = 1.0/_contactAssort.sumAllElements();
    _contactAssort = z * _contactAssort;
}




vector< vector<uint> > Simulator::draw_contacted_S_age_constraint(uint k,
                                                                  vector<uint> n_contacts,
                                                                  string infectious_type){
    /// Select the positions in '_indiv_S'
    /// for susceptibles that will be
    /// contacted by infectious individuals of a
    /// given infectious type (e.g. symptomatic or not).
    
    /// Selection is random, but with a constraint on ages.
    
    /// The returned format is a vector of vector, where
    /// column i can be seen as the list of susceptible
    /// contacted by the ith infector:
    ///
    /// infector_1     infector_2     infector_3    ...
    /// --------------------------------------------------
    /// indiv_S[.]     indiv_S[.]     indiv_S[.]    ...
    /// indiv_S[.]                    indiv_S[.]    ...
    /// indiv_S[.]                    indiv_S[.]    ...
    /// indiv_S[.]                                  ...
    
    
    uint n_inf = 0;
    vector<individual*> infector;
    vector<individual*> susceptible = _world[k]._indiv_S;
    
    uint n_susc = (uint)susceptible.size();
    
    if(infectious_type == "Is"){
        n_inf = (uint)_world[k]._indiv_Is.size();
        infector = _world[k]._indiv_Is;
    }
    else if(infectious_type == "Ia"){
        n_inf = (uint)_world[k]._indiv_Ia.size();
        infector = _world[k]._indiv_Ia;
    }
    
    // Some checks:
    stopif(sumElements(n_contacts)>n_susc, "Asking to draw more susceptibles than it exists!");
    stopif(n_contacts.size() != n_inf, "Number of contacts not consistent with number of infecious individuals.");
    
    vector< vector<uint> > selected_S(n_inf);
    
    double lambda = _modelParam.get_prm_double("contactAssort_lambda") ;
    std::exponential_distribution<double> expdist(lambda);

    
    for(uint i=0; i<n_inf; i++){
        
        uint age_inf = (uint)infector[i]->get_age();
        
        // Retrieve age of all susceptibles
        // and score them according to the
        // desired assortativity in '_contactAssort':
        vector<double> scores(n_susc);
        vector<double> age_s(n_susc); // DELETE when finish debug
        for(uint s=0; s<n_susc; s++)
        {
            uint age_susc = (uint)susceptible[s]->get_age();
            scores[s]     = 1.0 - _contactAssort(age_inf, age_susc);
            
            age_s[s] = susceptible[s]->get_age();
            
        }
        // Sort the scores (highest = better chance of contact)
        // but return the _original_ vector indexes:
        vector<size_t> u = sort_indexes(scores);
        
        // Pick randomly the rank (position in 'u').
        // The closer we want to be from the desired
        // contact pattern (given by _contactAssort)
        // the larger should the intensity lambda be:
        // (expected position = 1/lambda ; first position is prefered age)
        
        for(uint j=0; j<n_contacts[i]; j++)
        {
            uint rank = (uint)expdist(_RANDOM_GENERATOR);
            if (rank>= u.size()) rank = (uint)u.size()-1;
            
            stopif(rank >= n_susc, "Susceptible selected for assortative mixing does not exist.");
            
            selected_S[i].push_back((uint)u[rank]);
        }
    }
    return selected_S;
}




vector< vector<uint> > Simulator::transmission_attempts(uint k,
                                                         vector< vector<uint> > selected_S,
                                                         string infectious_type){
    /// Attempts transmission on all selected susceptible individuals
    
    uint n_ss = (uint) selected_S.size();
    vector< vector<uint> > transm_success(n_ss);
    
    bool homog = _modelParam.get_prm_bool("homogeneous_contact");
    
    if(homog){
        // Homogeneous contact: for testing purpose only. (but keep it!)
        // Always success because the contact rate
        // is understood as an _effective_ one in the homogeneous case:
        for (uint i=0; i<n_ss; i++)
            transm_success[i].resize(selected_S[i].size(), true);
    }
    
    if(!homog){
        std::uniform_real_distribution<float> unif01(0.0, 1.0);
        
        for (uint i=0; i<n_ss; i++)
        {
            transm_success[i].resize(selected_S[i].size(), false);
            
            for (uint j=0; j < selected_S[i].size(); j++){
                // TO DO: could be optimized:
                // because sampling with replacement,
                // we may attempt transmission on a susceptible
                // who had successful transmission previously
                
                double p_ij = -999.99;
                
                if (infectious_type == "Is")
                    p_ij = calc_proba_transmission(_world[k]._indiv_Is[i],
                                                   _world[k]._indiv_S[selected_S[i][j]]);
                if (infectious_type == "Ia")
                    p_ij = calc_proba_transmission(_world[k]._indiv_Ia[i],
                                                   _world[k]._indiv_S[selected_S[i][j]]);
                
                stopif(p_ij <0 && p_ij > 1.0, "Probability of transmissiom not b/w 0 and 1");
                
                if (unif01(_RANDOM_GENERATOR)<p_ij) {
                    transm_success[i][j] = true;
                }
            }
        }
    }
    return transm_success;
}


void Simulator::transmission_wiw(int k,
                                  vector<vector<uint> > selected_S,
                                  vector<vector<uint> > transm_success,
                                  string infectious_type){
    /// Records who infected who.
    
    for (uint i=0; i<selected_S.size(); i++) {
        for (uint s=0; s<selected_S[i].size(); s++) {
            // if transmission is successful:
            if (transm_success[i][s])
            {
                // Generation intervals
                
                _world[k]._indiv_S[selected_S[i][s]]->set_acquisition_time(_current_time);
                ID id_S = _world[k]._indiv_S[selected_S[i][s]]->get_id();
                float ti = __UNDEFINED_FLOAT;
                
                if(infectious_type=="Is"){
                    _world[k]._indiv_Is[i]->increment_num_secondary_cases();
                    _world[k]._indiv_Is[i]->push_ID_secondary_cases(id_S);
                    ti = _world[k]._indiv_Is[i]->get_acquisition_time();
                }
                
                if(infectious_type=="Ia"){
                    _world[k]._indiv_Ia[i]->increment_num_secondary_cases();
                    _world[k]._indiv_Ia[i]->push_ID_secondary_cases(id_S);
                    ti = _world[k]._indiv_Ia[i]->get_acquisition_time();
                }
                _world[k]._indiv_S[selected_S[i][s]]->set_acquisition_time_infector(ti);
                
                // Ages:
                
                if(infectious_type=="Is") _wiw_ages[0].push_back(_world[k]._indiv_Is[i]->get_age());
                if(infectious_type=="Ia") _wiw_ages[0].push_back(_world[k]._indiv_Ia[i]->get_age());
                _wiw_ages[1].push_back(_world[k]._indiv_S[selected_S[i][s]]->get_age());
                
                
            } // end-if-transmission-success
        } //end-for-s
    } //end-for-i
}

uint Simulator::transmission_activation(int k,
                                         vector< vector<uint> > selected_S,
                                         vector< vector<uint> > transm_success){
    /// Activate the succesful transmission attempts:
    
    uint incidence = 0;
    
    // Duplicates and not ordered 'selected_S'
    // messes up the deletion of newly infected.
    // Hence, need to clean-up this vector
    // before erasing elements of '_indiv_S':
    vector<uint> SS_melt_all = melt(selected_S);
    vector<uint> TS_melt = melt(transm_success);
    
    // Pick up only the selected susceptible
    // where the transmission was successfull:
    vector<uint> SS_melt;
    for (uint i=0; i<SS_melt_all.size(); i++) {
        if(TS_melt[i]) SS_melt.push_back(SS_melt_all[i]);
    }
    // clean-up: sort descending and remove duplicates
    // (important for easy deletions!)
    vector<uint> SS_melt_success = sort_remove_duplicate(SS_melt, false);
    unsigned long nss =SS_melt_success.size();
    
    // * * * WARNING * * *
    // for performance, put this here instead of inside for loop below,
    // but this assumes only ONE disease
    disease flu ;
    if(nss>0) flu = _world[k]._indiv_S[SS_melt_success[0]]->get_disease();
    
    std::uniform_real_distribution<float> unif_01(0.0, 1.0);
    
    // activate transmission:
    for (uint m=0; m<nss; m++) {
        
        // acquire the disease and (inside fct) set DOL & DOI
        _world[k]._indiv_S[SS_melt_success[m]]->acquireDisease();
        
        // determine if symptomatic:
        float imm_cell = _world[k]._indiv_S[SS_melt_success[m]]->get_immunity_cell();
        float frailty  = _world[k]._indiv_S[SS_melt_success[m]]->get_frailty();
        double p_sympt = calc_proba_symptomatic(imm_cell, frailty);
        
        bool is_symptomatic = ( unif_01(_RANDOM_GENERATOR) < p_sympt );
        if (is_symptomatic) {
            _world[k]._indiv_S[SS_melt_success[m]]->set_is_symptomatic(true);
            _world[k]._indiv_S[SS_melt_success[m]]->set_was_symptomatic(true);
            
            // For performance, the time and duration of hospitalization
            // is decided (in advance) at infection
            // (this is faster than checking if we hospitalized at every time step)
            double  p_hosp  = calc_proba_hospitalized(frailty);
            bool    willbe_hosp = ( unif_01(_RANDOM_GENERATOR) < p_hosp );
            if (willbe_hosp) {
                _world[k]._indiv_S[SS_melt_success[m]]->futureHospitalization();
            }
            
            // Also decide here if death at
            // the end of hospitalization period:
            double  p_death  = calc_proba_death(frailty);
            bool    will_die = ( unif_01(_RANDOM_GENERATOR) < p_death );
            if (will_die){
                _world[k]._indiv_S[SS_melt_success[m]]->futureDeath();
            }
        }
        
        // Book keeping: Update counters
        _world[k].update_epidemic_count(*_world[k]._indiv_S[SS_melt_success[m]], "new_case");
        
        // Book keeping: Update pointer tables:
        uint tmp1 = SS_melt_success[m];
        _world[k]._indiv_S.erase(_world[k]._indiv_S.begin() + tmp1);
        
        incidence ++;
        // TO DO: record infection times, etc.
    }
    return incidence;
}



uint Simulator::transmission_process(uint k, double dt, string infectious_type){
    
    stopif(k >= _world.size(), "Asking for an inexistent social place.");
    
    // If there are no susceptibles
    // or no infectious individuals
    // in the kth social place,
    // then no transmission can occur!!!
    uint nI = (infectious_type=="Is")?_n_Is:_n_Ia;
    if (_n_S == 0 || (nI==0) ) {
        return 0;
    }
    
    // Number of contacts for each infectious individual.
    // (numbers are stored in the order of vectors _indiv_Is, _indiv_Ia)
    vector<uint> n_contacts_Ix = draw_n_contacts(k, dt, infectious_type);
    
    // Randomly select susceptibles (position in '_indiv_S')
    // in this social place that will
    // be in contact: sampling _with_ replacement
    // (a S can be in contact with more than one I)
    
    vector< vector<uint> > selected_S_Ix = draw_contacted_S_age_constraint(k, n_contacts_Ix, infectious_type);
    //vector< vector<uint> > selected_S_Ix = draw_contacted_S(k, n_contacts_Ix, infectious_type);
    
    // Transmission attempts
    // (need to do this intermediary step
    // in order not to mess up pointers vector '_indiv_X')
    vector< vector<uint> > transm_success_Ix = transmission_attempts(k, selected_S_Ix, infectious_type);
    
    // Who infected who (and when)
    transmission_wiw(k, selected_S_Ix, transm_success_Ix, infectious_type);
    
    // Activate disease acquisition:
    uint inc_Ix = transmission_activation(k, selected_S_Ix, transm_success_Ix);
    
    return inc_Ix;
}


uint Simulator::transmission_oneSP(uint k,
                                    double dt){
    /// Performs transmission within the k^th social place.
    /// Returns incidence for _this_ social place, during the time step 'dt'.
    
    
    uint inc_Is, inc_Ia;
    
    // During a 'transmission_process'
    // susceptibles get infected and hence
    // are removed from the list of susceptibles '_indiv_S'.
    // So the order of running the process matters
    // (the firsts get more chance to infect as there are more
    // susceptible available).
    // Because we don't want to bias infections towards
    // either symptomatic or asymptomatic, the order is
    // drawn at random.
    
    std::uniform_real_distribution<> unif(0.0, 1.0);
    if (unif(_RANDOM_GENERATOR) < 0.5){
        inc_Is = transmission_process(k, dt, "Is");
        inc_Ia = transmission_process(k, dt, "Ia");
    }
    else{
        inc_Ia = transmission_process(k, dt, "Ia");
        inc_Is = transmission_process(k, dt, "Is");
    }
    return inc_Is + inc_Ia;
}


void Simulator::transmission_world(double timeslice){
    /// Simulates disease transmissions in the whole world (all social places)
    
    uint incidence = 0;
    for(uint k=0; k < _world.size(); k++){
        incidence += transmission_oneSP(k, timeslice);
    }
    _incidence = incidence;
}


uint Simulator::census_total_alive(){
    /// Counts all individuals that are alive
    uint cnt = 0;
    for(int k=0; k<_world.size(); k++) cnt += _world[k].census_alive();
    return cnt;
}


uint Simulator::prevalence(){
    /// Prevalence in the whole world
    
    uint cnt = 0;
    for (int k=0; k<_world.size(); k++) cnt += _world[k].get_prevalence();
    return cnt;
}


uint Simulator::population_size(){
    
    uint s = 0;
    for(int i=0; i<_world.size(); i++) s+=_world[i].get_size();
    return s;
}

dcDataFrame Simulator::census_sp(){
    
    vector<double> sp_id(_world.size());
    vector<double> sp_type(_world.size());
    vector<double> sp_nlinked(_world.size());
    
    
    for(uint k=0; k<_world.size(); k++){
        sp_id[k]      = (double)_world[k].get_id_sp();
        sp_type[k]    = (double)_world[k].get_type();
        sp_nlinked[k] = (double)_world[k].get_linked_indiv_id().size();
    }
    
    dcDataFrame x(sp_id,"sp_id");
    x.addcol("sp_type", sp_type);
    x.addcol("nlinked", sp_nlinked);
    return x;
}

void Simulator::display_summary_info(){
    
    
    cout << endl << endl;
    cout << " =======  SIMULATOR INFO ======= " << endl << endl;
    
    // Simulation parameters:
    
    tabcout("Simulation horizon (days)",_horizon);
    
    // Number of social places by type and
    // number of indiv linked to these SP:

    vector<uint> cnt(SP_MAX,0);
    vector<uint> cnt_indiv(SP_MAX,0);
    uint cnt_children = 0;
    
    for(uint k=0; k<_world.size(); k++){
        uint spt = _world[k].get_type();
        cnt[spt]++;
        cnt_indiv[spt] += (uint)_world[k].get_linked_indiv_id().size();
        cnt_children   += _world[k].census_alive_age(0, 17.99999999
                                                     );
    }
    
    cout << endl << "Number of social places by type:" <<endl;
    for (uint i=0; i<SP_MAX; i++){
        tabcout(SPtype2string((SPtype)i), cnt[i]);
    }
    stopif(sumElements(cnt) != _world.size(),
           "Total number of social places not consistent.");
    tabcout("TOTAL", _world.size());
    
    // Number of individuals:
    cout << endl << "Number of individual linked to social place type:" <<endl;
    for (uint i=0; i<SP_MAX; i++){
        tabcout(SPtype2string((SPtype)i), cnt_indiv[i]);
    }
    cout << endl;
    tabcout("Total number of individuals", world_size(_world));
    tabcout("Total number of children", cnt_children);
    
    // Schedules used:
    cout << endl;
    tabcout("Individuals with schedule worker_sed", census_schedule(_world, "worker_sed"), 40);
    tabcout("Individuals with schedule student",    census_schedule(_world, "student"), 40);
    tabcout("Individuals with schedule unemployed", census_schedule(_world, "unemployed"), 40);
    
    // Interventions:
    cout << endl << "Number of interventions : " << _intervention.size() <<endl;
    for(int i=0; i<_intervention.size(); i++) {
        cout << " Intervention["<<i<<"]";
     _intervention[i].display_info();
    }
    
    cout << endl;
    cout << " = (end of simulator info) = " << endl << endl;
}


void Simulator::display_split_pop_present(){
    
    cout<<endl<<"------------"<<endl;
    uint s = 0;
    uint p = 0;
    for(int i=0;i<_world.size();i++){
        cout << "sp_"<<i<<" : present = " << _world[i].get_size()<<" (prev="<<_world[i].get_prevalence()<<")";
        cout << "\t ["<< SPtype2string(_world[i].get_type())<<"]" <<endl;
        s += _world[i].get_size();
        p += _world[i].get_prevalence();
    }
    cout<< "Total population = "<< s << " (prev="<<p<<")"<<endl;
    cout<<"------------"<<endl;
    
    
}

void Simulator::display_split_pop_linked(){
    cout<<endl<<"------------"<<endl;
    uint s = 0;
    for(int i=0;i<_world.size();i++){
        cout << "sp_"<<i<<" : Linked indiv = " << _world[i].n_linked_indiv();
        cout << "\t ["<< SPtype2string(_world[i].get_type())<<"]" <<endl;
        s += _world[i].n_linked_indiv();
    }
    cout<< "Total linked = "<< s << endl;
    cout<<"------------"<<endl;
}


void Simulator::seed_infection(vector<ID> id_sp, vector<uint> I0){
    
    stopif(id_sp.size() != I0.size(), "vectors must be same size");
    
    vector<vector<uint> > selected_S(1);
    vector<vector<uint> > transm_success(id_sp.size());
    
    for(ID cnt=0; cnt<id_sp.size();cnt++)
    {
        for(ID i=0; i<_world.size(); i++)
        {
            if (_world[i].get_id_sp() == id_sp[cnt])
            {
                cout << " Seeding infection in SP_" << id_sp[cnt] << endl;
                
                // Consistency check:
                string errmsg =  "Cannot seed in SP_ID_" + to_string(i) + " because it has less individuals than intial infections requested!";
                stopif(_world[i].get_size() < I0[cnt], errmsg);
                
                // Infection:
                for(uint j=0; j<I0[cnt]; j++) {
                    selected_S[0].push_back(j);
                    transm_success[0].push_back(true);
                }
                transmission_activation(id_sp[cnt], selected_S, transm_success);
            }
        }
    }
    update_pop_count();
}


void Simulator::displayInfo_indiv(){
    /// Display informations on all individuals, in all social places:
    
    for (uint k=0; k<_world.size(); k++) {
        for (ID i=0; i<_world[k].get_size(); i++) {
            _world[k].get_indiv(k).displayInfo();
        }
    }
}


dcDataFrame Simulator::timeseries(){
    
    dcDataFrame df(_ts_times,"time");
    
    df.addcol("incidence", to_vector_double(_ts_incidence));
    df.addcol("prevalence", to_vector_double(_ts_prevalence));
    df.addcol("nS", to_vector_double(_ts_S));
    df.addcol("nE", to_vector_double(_ts_E));
    df.addcol("nIa", to_vector_double(_ts_Ia));
    df.addcol("nIs", to_vector_double(_ts_Is));
    df.addcol("nR", to_vector_double(_ts_R));
    df.addcol("nD", to_vector_double(_ts_D));
    df.addcol("n_treated", to_vector_double(_ts_n_treated));
    df.addcol("n_vaccinated", to_vector_double(_ts_n_vaccinated));
    
    return df;
}



void Simulator::update_ts_census_SP(){
    /// Update the data frame recording
    /// the number of individuals in each
    /// disease stage, for every social places.
    
    unsigned long n = _world.size();
    
    for (uint k=0; k<n; k++) {
        _ts_census_sp_time.push_back(_current_time);
        _ts_census_sp_id.push_back(_world[k].get_id_sp());
        _ts_census_sp_type.push_back( SPtype2string(_world[k].get_type()));
        _ts_census_sp_nS.push_back(_world[k].get_n_S());
        _ts_census_sp_nE.push_back(_world[k].get_n_E());
    }
}



void Simulator::update_pop_count(){
    /// Update the count of individuals in a all
    /// stages of the disease at the 'simulation' level
    /// (social places were updated during transmission).
    
    unsigned long n = _world.size();
    _n_S  = 0;
    _n_E  = 0;
    _n_Ia = 0;
    _n_Is = 0;
    _n_R  = 0;
    
    for (ID i=0; i<n; i++) {
        _n_S  += _world[i].get_n_S();
        _n_E  += _world[i].get_n_E();
        _n_Ia += _world[i].get_n_Ia();
        _n_Is += _world[i].get_n_Is();
        _n_R  += _world[i].get_n_R();
    }
    _prevalence = _n_E + _n_Ia + _n_Is;
}


void Simulator::check_book_keeping(){
    /// Check consistency of book keeping.
    /// WARNING: FOR DEBUG ONLY, SLOWS DOWN EXECUTION!

    check_sp_integrity(_world);

    // more individual present than linked?
    for (ID k=0; k<_world.size(); k++) {
        unsigned long pres = _world[k].get_indiv().size();
        unsigned long linked = _world[k].get_linked_indiv_id().size();
        stopif(pres>linked, "More individuals present than linked in SP ID "+to_string(k));
    }
    
    // stage S
    uint nS=0;
    uint nS_census=0;
    for (ID k=0; k<_world.size(); k++) {
        nS += _world[k].get_n_S();
        nS_census += _world[k].census_disease_stage("S");
        stopif(_world[k].get_n_S() != _world[k].get_id_S().size() ||
               _world[k].get_n_S() != _world[k]._indiv_S.size(),
               "Book keeping error with S IDs.");
        
        for(uint i=0; i<_world[k]._indiv_S.size(); i++){
            ID id = _world[k]._indiv_S[i]->get_id();
            uint j=0;
            while(_world[k].get_indiv(j).get_id() != id)
                j++;
            stopif(j >= _world[k].get_indiv().size(), "Book keeping error with _indiv_S.");
            stopif(!_world[k].get_indiv(j).is_susceptible(), "Book keeping error with _indiv_S.");
        }
    }
    bool check_S = ( (_n_S == nS) && (nS == nS_census) );
    stopif(!check_S, "Book keeping error with S stage");
    
    // stage E
    uint nE=0;
    uint nE_census=0;
    for (ID k=0; k<_world.size(); k++) {
        nE += _world[k].get_n_E();
        nE_census += _world[k].census_disease_stage("E");
    }
    bool check_E = ( (_n_E == nE) && (nE == nE_census) );
    stopif(!check_E, "Book keeping error with E stage");
    
    // stage Is
    uint nIs=0;
    uint nIs_census=0;
    for (ID k=0; k<_world.size(); k++) {
        nIs += _world[k].get_n_Is();
        nIs_census += _world[k].census_disease_stage("Is");
        stopif(_world[k].get_n_Is() != _world[k].get_id_Is().size() ||
               _world[k].get_n_Is() != _world[k]._indiv_Is.size(),
               "Book keeping error with Is IDs.");
        
        for(uint i=0; i<_world[k]._indiv_Is.size(); i++){
            ID id = _world[k]._indiv_Is[i]->get_id();
            uint j=0;
            while(_world[k].get_indiv(j).get_id() != id)
                j++;
            stopif(j >= _world[k].get_indiv().size(), "Book keeping error with _indiv_Is.");
            stopif(!_world[k].get_indiv(j).is_infectious(), "Book keeping error with _indiv_Is.");
            stopif(!_world[k].get_indiv(j).is_symptomatic(), "Book keeping error with _indiv_Is.");
        }
        
    }
    bool check_Is = ( (_n_Is == nIs) && (nIs == nIs_census) );
    stopif(!check_Is, "Book keeping error with Is stage");
    
    // stage Ia
    uint nIa=0;
    uint nIa_census=0;
    for (ID k=0; k<_world.size(); k++) {
        nIa += _world[k].get_n_Ia();
        nIa_census += _world[k].census_disease_stage("Ia");
        stopif(_world[k].get_n_Ia() != _world[k].get_id_Ia().size() ||
               _world[k].get_n_Ia() != _world[k]._indiv_Ia.size(),
               "Book keeping error with Ia IDs.");
        
        for(uint i=0; i<_world[k]._indiv_Ia.size(); i++){
            ID id = _world[k]._indiv_Ia[i]->get_id();
            uint j=0;
            while(_world[k].get_indiv(j).get_id() != id)
                j++;
            stopif(j >= _world[k].get_indiv().size(), "Book keeping error with _indiv_Ia.");
            stopif(!_world[k].get_indiv(j).is_infectious(), "Book keeping error with _indiv_Ia.");
            stopif(_world[k].get_indiv(j).is_symptomatic(), "Book keeping error with _indiv_Ia.");
        }
        
    }
    bool check_Ia = ( (_n_Ia == nIa) && (nIa == nIa_census) );
    stopif(!check_Ia, "Book keeping error with Ia stage");
    
    // stage R
    uint nR=0;
    uint nR_census=0;
    for (ID k=0; k<_world.size(); k++) {
        nR += _world[k].get_n_R();
        nR_census += _world[k].census_disease_stage("R");
    }
    bool check_R = ( (_n_R == nR) && (nR == nR_census) );
    stopif(!check_R, "Book keeping error with R stage");
    
    // stage H
    uint nH=0;
    uint nH_census=0;
    for (ID k=0; k<_world.size(); k++) {
        nH += _world[k].get_n_H();
        nH_census += _world[k].census_disease_stage("H");
        
        for(uint i=0; i<_world[k]._indiv_H.size(); i++){
            ID id = _world[k]._indiv_H[i]->get_id();
            uint j=0;
            while(_world[k].get_indiv(j).get_id() != id)
                j++;
            stopif(j >= _world[k].get_indiv().size(), "Book keeping error with _indiv_H.");
            stopif(!_world[k].get_indiv(j).is_hosp(), "Book keeping error with _indiv_H.");
        }
    }
    bool check_H = ( (_n_H == nH) && (nH == nH_census) );
    stopif(!check_H, "Book keeping error with H stage");
}



void Simulator::define_all_id_tables(){
    
    for (uint k=0; k<_world.size(); k++)
    {
        vector<individual> indiv = _world[k].get_indiv();
        
        _world[k].clear_id_S();
        _world[k].clear_id_Is();
        _world[k].clear_id_Ia();
        _world[k]._indiv_S.clear();
        _world[k]._indiv_Is.clear();
        _world[k]._indiv_Ia.clear();
        _world[k]._indiv_H.clear();
        _world[k]._indiv_vax.clear();
        
        for (uint i=0; i< indiv.size(); i++)
        {
            if (indiv[i].is_susceptible()) {
                _world[k].add_id_S(indiv[i].get_id());
                _world[k]._indiv_S.push_back(_world[k].get_mem_indiv(i));
            }
            else {
                if (indiv[i].is_infectious() && indiv[i].is_symptomatic()){
                    _world[k].add_id_Is(indiv[i].get_id());
                    _world[k]._indiv_Is.push_back(_world[k].get_mem_indiv(i));
                }
                else if (indiv[i].is_infectious() && !indiv[i].is_symptomatic()){
                    _world[k].add_id_Ia(indiv[i].get_id());
                    _world[k]._indiv_Ia.push_back(_world[k].get_mem_indiv(i));
                }
                if (indiv[i].is_hosp() && indiv[i].is_alive()){
                    _world[k]._indiv_H.push_back(_world[k].get_mem_indiv(i));
                }
            }
            if  (indiv[i].is_vaccinated()){
                _world[k]._indiv_vax.push_back(_world[k].get_mem_indiv(i));
            }
        }
    }
}


double Simulator::select_contact_rate_ratio(double age,
                                            SPtype sp_type){
    
    double ratio_indiv = 1.0;
    double ratio_sp    = 1.0;
    
    // age:
    if (1 <   age && age < 10) ratio_indiv = ratio_indiv * _modelParam.get_prm_double("contact_ratio_age_1_10");
    if (10 <= age && age < 16) ratio_indiv = ratio_indiv * _modelParam.get_prm_double("contact_ratio_age_10_16");
    if (65 <= age)             ratio_indiv = ratio_indiv * _modelParam.get_prm_double("contact_ratio_age_over_65");
    
    // social place:
    if (sp_type == SP_household) ratio_sp = ratio_sp * _modelParam.get_prm_double("contact_ratio_sp_household");
    if (sp_type == SP_pubTransp) ratio_sp = ratio_sp * _modelParam.get_prm_double("contact_ratio_sp_pubTransport");
    if (sp_type == SP_school)    ratio_sp = ratio_sp * _modelParam.get_prm_double("contact_ratio_sp_school");
    
    return ratio_indiv * ratio_sp;
}


double  Simulator::draw_contact_rate(individual* indiv, uint k){

    
    bool homog = _modelParam.get_prm_bool("homogeneous_contact");
    
    double cr = __UNDEFINED_DOUBLE;
    
    if(homog){
        // Homogeneous contacts.
        // Used for testing purposes (e.g. compare to ODE models).
        uint ns = _world[k].get_n_S();
        uint n  = (uint)_world[k].get_size();
        cr      = _modelParam.get_prm_double("contact_rate_mean") * ns / n;
    }
    
    if(!homog){
        // Contact rate depends on features of
        // both individual and social place
        SPtype sp_type = _world[k].get_type();
        double age     = indiv->get_age();
        
        double ratio = select_contact_rate_ratio(age, sp_type);
        
        // The contact rate is drawn from another distribution
        // in order to allow for super-spreading events
        double cr_mean = ratio * _modelParam.get_prm_double("contact_rate_mean");
        double cr_sd   = _modelParam.get_prm_double("contact_rate_stddev");
        
        //DELETE: std::exponential_distribution<double> expdist(1.0 / (ratio * cr_baseline));
        double tmp  = 1 + cr_sd*cr_sd/cr_mean/cr_mean;
        double m_ln = log(cr_mean/sqrt(tmp));
        double s_ln = sqrt(log(tmp));
        std::lognormal_distribution<double> lndist(m_ln, s_ln);
        
        cr = lndist(_RANDOM_GENERATOR);
        
        
        //DEBUG:
        //cout<<"cr_mean= " << cr_mean << " ; ratio = "<< ratio << " ; cr = " << cr << endl;
    }
    return cr;
}



void Simulator::assign_hospital_to_individuals(){
    /// Assign hospitals to individuals
    
    // TO DO: implement something more elaborated
    // (based on distance, when coordinate available)
    
    // SIMPLE ALGO FOR TEST = assign everyone to the first hospital
    
    //find first hospital
    uint k_hosp = 0;
    while (_world[k_hosp].get_type() != SP_hospital && k_hosp<_world.size()){
        k_hosp++;
    }
    stopif(k_hosp>=_world.size(), "No hospital found!");
    
    for (uint k=0; k<_world.size(); k++) {
        for (uint i=0; i<_world[k].get_size() ; i++) {
            _world[k].set_id_sp_hospital(i, _world[k_hosp]);
        }
    }
}


void Simulator::hospitalize(){
    
    std::uniform_real_distribution<> unif(0.0, 1.0);
    
    for (uint k=0; k<_world.size(); k++)
    {
        if(_world[k].get_type() != SP_hospital)
        {
            uint nIs = _world[k].get_n_Is();
            
            stopif(nIs != _world[k]._indiv_Is.size(), "Book keeping issue with Is!");
            
            // * * * IMPORTANT * * *
            // Must run this loop
            // in _descending_ order, else
            // it messes up the pointers vector deletion
            for (uint i=nIs;  (i--) >0;) {
                
                // If time to hospitalize:
                if (!_world[k]._indiv_Is[i]-> is_hosp() &&
                    _world[k]._indiv_Is[i]-> willbe_hosp() &&
                    _world[k]._indiv_Is[i]-> time_to_hospitalize() &&
                    !_world[k]._indiv_Is[i]->is_discharged())
                {
                    // set the hospitalization flag for this individual
                    _world[k]._indiv_Is[i]->set_is_hosp(true);
                    _world[k]._indiv_Is[i]->set_doh(SUPERTINY);
                    
                    ID id_sp_hospital = _world[k]._indiv_Is[i]->get_id_sp_hospital();
                    ID id_indiv = _world[k]._indiv_Is[i]->get_id();
                    uint pos_indiv = _world[k].find_indiv_pos(id_indiv);
                    
                    stopif(id_sp_hospital==__UNDEFINED_ID,
                           "Individual ID "+to_string(id_indiv)+
                           " is not linked to a hospital. " +
                           "Consider increasing the size of the hospital(s).");
                    
                    // update counters
                    _n_H++;
                    _world[id_sp_hospital].increment_n_H();
                    
                    // Move individual to its linked SP_hospital
                    move_one_individual(pos_indiv, k, id_sp_hospital);
                    // cout <<" DEBUG: indiv ID_"<<to_string(id_indiv)<<"  hospitalized" <<endl;
                }
            }
        }
    }
}


void Simulator::discharge_hospital(uint idx_timeslice){
    
    for (uint k=0; k<_world.size(); k++)
    {
        uint nH = _world[k].get_n_H();
        
        stopif(nH != _world[k]._indiv_H.size(), "Book keeping issue with H!");
        
        if(_world[k].get_type() == SP_hospital &&   // <-- can only leave from a hospital
           _world[k]._indiv_H.size() >0 )           // <-- at least one must be hospitalized
        {
            // * * * IMPORTANT * * *
            // Must run this loop
            // in _descending_ order, else
            // it messes up the pointers vector deletion
            for (uint i=nH;  (i--) >0;) {
                
                // If time to leave hospital:
                if ( _world[k]._indiv_H[i]-> is_discharged() &&
                     _world[k]._indiv_H[i]-> is_alive()      &&
                    !_world[k]._indiv_H[i]-> will_die()
                    )
                {
                    
                    ID id_sp_hospital = _world[k]._indiv_H[i]->get_id_sp_hospital();
                    ID id_indiv       = _world[k]._indiv_H[i]->get_id();
                    uint pos_indiv    = _world[k].find_indiv_pos(id_indiv);
                    
                    // set the hospitalization flag for this individual
                    _world[k]._indiv_H[i]->set_is_hosp(false);
                    _world[k]._indiv_H[i]->set_was_hosp(true);
                    
                    // update counters
                    _n_H--;
                    _world[id_sp_hospital].decrement_n_H();
                    
                    // Move individual to its scheduled social place.
                    // Retrieve its actual destination:
                    ID id_dest = _world[k].find_dest(i, idx_timeslice);
                    
                    stopif(id_dest==__UNDEFINED_ID,
                           "Discharge hopital: Undefined destination for indiv ID_" + to_string(id_indiv));
                    
                    if ( id_dest!=k )
                    { move_one_individual(pos_indiv, k, id_dest); }
                    
                    //                    cout <<" DEBUG: indiv ID_"<<to_string(id_indiv)<<" leave hospital!" <<endl;
                }
            }
        } // end-if-SP_hospital
    }
}


void Simulator::death_hospital(){
    
    /// Scan all infectious hospitalized individuals
    /// and trigger death if meant to die (see '_will_die')
    
    for (uint k=0; k<_world.size(); k++)
    {
        uint nH = _world[k].get_n_H();
        
        stopif(nH != _world[k]._indiv_H.size(), "Book keeping issue with H!");
        
        if(_world[k].get_type() == SP_hospital &&   // <-- can only leave from a hospital
           _world[k]._indiv_H.size() >0 )           // <-- at least one must be hospitalized
        {
            // * * * IMPORTANT * * *
            // Must run this loop
            // in _descending_ order, else
            // it messes up the pointers vector deletion
            for (uint i=nH;  (i--) >0;) {
                
                // If time to die:
                if (!_world[k]._indiv_H[i]-> is_discharged() &&
                    _world[k]._indiv_H[i]-> is_alive()      &&
                    _world[k]._indiv_H[i]-> will_die()
                    )
                {
                    ID id_sp_hospital = _world[k]._indiv_H[i]->get_id_sp_hospital();
                    ID id_indiv       = _world[k]._indiv_H[i]->get_id();
                    //                    uint pos_indiv    = _world[k].find_indiv_pos(id_indiv);
                    
                    // set the hospitalization flag for this individual
                    _world[k]._indiv_H[i]->die();
                    
                    // update counters
                    _n_H--;
                    _n_Is--;
                    _n_D++;
                    _world[id_sp_hospital].decrement_n_H();
                    _world[id_sp_hospital].decrement_n_Is();
                    _world[id_sp_hospital].increment_n_D();
                    
                    // An hospitalized individual
                    // is in _both_ Is and H, so must
                    // remove it from these 2 categories:
                    removeValue(_world[id_sp_hospital]._indiv_Is, _world[k]._indiv_H[i]);
                    removeValue(_world[id_sp_hospital]._indiv_H,  _world[k]._indiv_H[i]); // <-- remove in _indiv_H after all other _indiv_X !!!
                    _world[id_sp_hospital].remove_id_Is(id_indiv);
                }
            }
        } // end-if-SP_hospital
    }
}



vector<individual*> Simulator::draw_targeted_individuals(uint i,
                                                         ID id_sp,
                                                         double dt){
    vector<individual*> indiv_drawn;
    bool found = false;
    
    float cvg_rate         = _intervention[i].get_cvg_rate();
    string type_target     = _intervention[i].get_type_indiv_targeted();
    
    if(type_target == "symptomatic"){
        
        // * * WARNING * *  excludes symptomatic already treated!
        
        found = true;
        uint nI = (uint)_world[id_sp]._indiv_Is.size();
        
        if (nI > 0){
            // Draw the total number of individuals that are targeted
            std::poisson_distribution<> poiss(cvg_rate * nI * dt);
            uint n_target = poiss(_RANDOM_GENERATOR);
            if (n_target > nI) n_target = nI;
            // Select targeted individuals:
            uint cnt = 0;
            for(uint i=0; i<nI && cnt<n_target; i++) {
                bool is_treated = _world[id_sp]._indiv_Is[i]->is_treated();
                if (!is_treated){
                    indiv_drawn.push_back(_world[id_sp]._indiv_Is[i]);
                    cnt ++;
                }
            }
        }
    }
    
    else if(type_target == "susceptible"){
        
        // * * WARNING * * excludes susceptible already treated & vaccinated!
        
        found = true;
        uint nS  = (uint)_world[id_sp]._indiv_S.size();
        
        if (nS > 0){
            // Draw the total number of individuals that are targeted
            
            float intensity = cvg_rate * nS * dt;
            std::poisson_distribution<> poiss(intensity);
            uint n_target = poiss(_RANDOM_GENERATOR);
            if (n_target > nS) n_target = nS;
            // Select targeted individuals:
            uint cnt = 0;
            for(uint j=0; j<nS && cnt<n_target; j++) {
                bool is_treated = _world[id_sp]._indiv_S[j]->is_treated();
                bool is_vax     = _world[id_sp]._indiv_S[j]->is_vaccinated();
                if (!is_treated && !is_vax){
                    indiv_drawn.push_back(_world[id_sp]._indiv_S[j]);
                    cnt ++;
                }
            }
        }
    } // end-if-type_target == "susceptible"
    
    else if(type_target == "young_old"){
        
        // * * WARNING * * targets _ALL_ (including infected) young/old individuals
        // as long as not already vaccinated or treated (if treated,
        // it is necessarily symptomatic, but not all symptomatic are treated)
    
        found = true;
        float age_old   = AGE_OLD;
        float age_young = AGE_YOUNG;
        
        // First, count how many individuals
        // of the targeted ages are present
        // in this social place:
        uint n_yo = 0;
        vector<uint> pos_yo;
        
        for (uint j=0; j< _world[id_sp].get_size(); j++) {
            double age = _world[id_sp].get_indiv(j).get_age();
            if(age <= age_young || age >= age_old) {
                n_yo++;
                pos_yo.push_back(j);
            }
        }
        
        if(n_yo>0){
            // Draw the total number of individuals
            // that are targeted:
            float intensity = cvg_rate * n_yo * dt;
            
            //DEBUG
//            cout << " DEBUG YO: SP_id = "<< id_sp << endl;
//            cout << " DEBUG YO: n_yo = "<< n_yo<< endl;
//            cout << " DEBUG YO: intensity = "<< intensity << endl;
            
            std::poisson_distribution<> poiss(intensity);
            uint n_target = poiss(_RANDOM_GENERATOR);
            if (n_target > n_yo) n_target = n_yo;
            
            // Select targeted individuals:
            uint cnt = 0;
            for(uint j=0; j<pos_yo.size() && cnt<n_target ; j++) {
                bool is_treated = _world[id_sp].get_indiv(pos_yo[j]).is_treated();
                bool is_vax     = _world[id_sp].get_indiv(pos_yo[j]).is_vaccinated();
                if (!is_treated && !is_vax){
                    individual* tmp = _world[id_sp].get_mem_indiv(pos_yo[j]);
                    indiv_drawn.push_back(tmp);
                    cnt ++;
                }
                //DEBUG
                //cout << "DEBUG:: YO vax = "<< cnt <<endl;
            }
        }
    }// end-if-type_target == "young_old"
    
    stopif(!found, "Type of targeted individual unknown: " + type_target);
    return indiv_drawn;
}


void Simulator::count_targeted_by_intervention(){
    
    _max_cvg_interv.resize(_intervention.size());
    
    for (uint i=0; i<_intervention.size(); i++){
        
        string type_target     = _intervention[i].get_type_indiv_targeted();
        
        uint cnt = 0;
        
        if(type_target == "symptomatic"){
            
            // WARNING : in this case, is is not possible
            // to know in advance the targeted population,
            // so max coverage is relative to total size (=susceptible at start)
            for(uint id_sp=0; id_sp< _world.size(); id_sp++)
                cnt += (uint)_world[id_sp]._indiv_S.size();
        }
        
        if(type_target == "susceptible"){
            for(uint id_sp=0; id_sp< _world.size(); id_sp++)
                cnt += (uint)_world[id_sp]._indiv_S.size();
        }
        
        if(type_target == "young_old"){
            float age_old   = AGE_OLD;
            float age_young = AGE_YOUNG;
            
            for(uint id_sp=0; id_sp< _world.size(); id_sp++){
                uint n_yo = 0;
                // counts individuals of the targeted age:
                for (uint i=0; i< _world[id_sp].get_size(); i++) {
                    double age = _world[id_sp].get_indiv(i).get_age();
                    if(age <= age_young || age >= age_old) {
                        n_yo++;
                    }
                }
                cnt += n_yo;
            }
            
        }
        // Maximum number of individuals
        // targeted by ith intervention:
        _max_cvg_interv[i] = (uint)(cnt * _intervention[i].get_cvg_max_proportion());
    }
}
    



void Simulator::activate_interventions(ID id_sp, double dt,
                                       float treat_doi_reduc,
                                       float vax_imm_hum_incr,
                                       float vax_imm_cell_incr,
                                       float vax_frail_incr,
                                       float vax_lag)
{
    _n_treated.resize(_intervention.size());
    _n_vaccinated.resize(_intervention.size());
    
    for (uint i=0; i<_intervention.size(); i++)
    {
        string interv_type = _intervention[i].get_type_intervention();
        
        if (_intervention[i].get_time_start() <= _current_time &&
            _intervention[i].get_time_end() > _current_time)
        {
            bool cond1 = (_n_treated[i]    > _max_cvg_interv[i]) && (interv_type=="treatment");
            bool cond2 = (_n_vaccinated[i] > _max_cvg_interv[i]) && (interv_type=="vaccination");
            
            if ( ! (cond1 || cond2) ){  // <-- check if max number targeted reached.
                
                vector<individual*> x = draw_targeted_individuals(i, id_sp, dt);
                
                // DEBUG:
//                if(x.size()>0){
//                    cout << " DEBUG: drawn size = "<< x.size()<< " from SP_id:"<< id_sp <<endl;
//                    
//                }
                
                _intervention[i].act_on_individual(x,
                                                   _current_time,
                                                   treat_doi_reduc,
                                                   vax_imm_hum_incr,
                                                   vax_imm_cell_incr,
                                                   vax_frail_incr,
                                                   vax_lag);
                
                if      (_intervention[i].get_type_intervention()=="treatment")
                    _n_treated[i] += x.size();
                
                else if (_intervention[i].get_type_intervention()=="vaccination")
                {
                    _n_vaccinated[i] += x.size();
                    _world[id_sp]._indiv_vax.insert(_world[id_sp]._indiv_vax.end(),
                                                    x.begin(),
                                                    x.end());
                }
            } // end-if con1 cond2
            else{
//                cout << "DEBUG: max targeted for intervention #"<<i<<" reached."<<endl;
//                cout << "_n_treat="<<_n_treated[i]   <<">"<<_max_cvg_interv[i]<<endl;
//                cout << "_n_vax="  <<_n_vaccinated[i]<<">"<<_max_cvg_interv[i]<<endl;
            }
        } // end-if time
    } // end-for intervention  size
}


void Simulator::update_immunities() {

    for (uint k=0; k<_world.size(); k++) {
        unsigned long n_vax = _world[k]._indiv_vax.size();
        for (uint i=0; i<n_vax; i++)
        {
            float t_vax   = _world[k]._indiv_vax[i]->get_vax_time_received();
            float lag_vax = _world[k]._indiv_vax[i]->get_vax_lag_full_efficacy();
            
            if (_current_time < t_vax+lag_vax){
                
                float imm_hum0   =  _world[k]._indiv_vax[i]->get_imm_hum_when_recv_vax();
                float imm_cell0  =  _world[k]._indiv_vax[i]->get_imm_cell_when_recv_vax();
                
                float target_imm_hum  =  _world[k]._indiv_vax[i]->get_vax_target_immunity_hum();
                float target_imm_cell =  _world[k]._indiv_vax[i]->get_vax_target_immunity_cell();

                // Linear growth to target level:
                float curr_imm_hum   = linear_interpol(_current_time,
                                                       t_vax, imm_hum0,
                                                       t_vax+lag_vax, target_imm_hum);
                float curr_imm_cell  = linear_interpol(_current_time,
                                                       t_vax, imm_cell0,
                                                       t_vax+lag_vax, target_imm_cell);
                // update levels:
                _world[k]._indiv_vax[i]->set_immunity_hum(curr_imm_hum);
                _world[k]._indiv_vax[i]->set_immunity_cell(curr_imm_cell);
            }
        }
    }
}



bool Simulator::at_least_one_vaccination_intervention(){
    /// Check if there is at least one vaccination among all
    /// interventions defined for this Simulator
    
    for (uint i=0; i<_intervention.size(); i++) {
        if (_intervention[i].get_type_intervention() == "vaccination")
            return true;
    }
    return false;
}



bool Simulator::at_least_one_infected(){
    /// Check if there is at least one infected individual
    /// in the whole population.
    bool res = false;
    for (uint k=0; k<_world.size(); k++) {
        if (_world[k].get_n_E() + _world[k].get_n_Ia() + _world[k].get_n_Is() >0) {
            return true;
        }
    }
    return res;
}




void Simulator::assign_dox_distribution(string dol_distrib,
                                        string doi_distrib,
                                        string doh_distrib){
    /// Assign distributions for duration of latency,
    /// infectiousness and hospitalization
    for (uint k=0; k< _world.size(); k++)
    {
        unsigned long nk = _world[k].get_size();
        for (uint i=0; i< nk; i++) {
            _world[k].set_dol_distrib(i, dol_distrib);
            _world[k].set_doi_distrib(i, doi_distrib);
            _world[k].set_doh_distrib(i, doh_distrib);
        }
    }
}


void Simulator::assign_immunity_hum(){
    
    double baseline = _modelParam.get_prm_double("imm_hum_baseline");
    double agezero  = _modelParam.get_prm_double("imm_hum_agezero");
    double p        = _modelParam.get_prm_double("imm_hum_p");
    
    for (uint k=0; k< _world.size(); k++)
    {
        unsigned long nk = _world[k].get_size();
        for (uint i=0; i< nk; i++) {
            float age = _world[k].get_indiv(i).get_age();
            double immh = immunity_humoral(age, agezero, baseline, p);
            _world[k].set_immunity_hum(i, immh);
        }
    }
}


void Simulator::assign_immunity_cell(){
    
    double imm_max = _modelParam.get_prm_double("imm_cell_max");
    double slope   = _modelParam.get_prm_double("imm_cell_slope");
    double pivot   = _modelParam.get_prm_double("imm_cell_pivot");
    
    for (uint k=0; k< _world.size(); k++)
    {
        unsigned long nk = _world[k].get_size();
        std::uniform_real_distribution<float> unif01(0,1);
        
        for (uint i=0; i< nk; i++) {
            float age = _world[k].get_indiv(i).get_age();
            double immc = immunity_cellular(age, imm_max, slope, pivot);
            _world[k].set_immunity_cell(i, immc);
        }
    }
}


void Simulator::assign_frailty(){
    /// Calculate frailty index for all individuals
    
    // Retrieve model parameters for frailty:
    float f0         = _modelParam.get_prm_double("frailty_0");
    float fmin       = _modelParam.get_prm_double("frailty_min");
    float agemin     = _modelParam.get_prm_double("frailty_agemin");
    float agepivot   = _modelParam.get_prm_double("frailty_agepivot");
    float fpivot     = _modelParam.get_prm_double("frailty_pivot");
    float powerChild = _modelParam.get_prm_double("frailty_powerChild");
    float frail_sd   = _modelParam.get_prm_double("frailty_sd");
    
    // loop through all indviduals:
    for (uint k=0; k< _world.size(); k++)
    {
        unsigned long nk = _world[k].get_size();
        for (uint i=0; i< nk; i++) {
            float age = _world[k].get_indiv(i).get_age();
            float frail_mean = frailty_mean(age, f0, fmin, agemin, agepivot, fpivot, powerChild);
            std::normal_distribution<> norm(frail_mean,frail_sd);
            
            float frail = norm(_RANDOM_GENERATOR);
            if(frail<0) frail = 0.0;
            else if(frail>1) frail = 1.0;
            
            _world[k].set_frailty(i, frail);
        }
    }
}



void Simulator::timeseries_update(){
    _ts_times.push_back(_current_time);
    _ts_incidence.push_back(_incidence);
    _ts_prevalence.push_back(_prevalence);
    _ts_S.push_back(_n_S);
    _ts_E.push_back(_n_E);
    _ts_Ia.push_back(_n_Ia);
    _ts_Is.push_back(_n_Is);
    _ts_H.push_back(_n_H);
    _ts_R.push_back(_n_R);
    _ts_D.push_back(_n_D);
    _ts_n_treated.push_back(sumElements(_n_treated));
    _ts_n_vaccinated.push_back(sumElements(_n_vaccinated));
}


void Simulator::change_rnd_sp_other(){    
    // Probability to change 'other' social place
    // during a time slice:
    double p = _modelParam.get_prm_double("proba_change_sp_other");
    
    for (uint k=0; k<_world.size(); k++) {
        unsigned long n_indiv = _world[k].get_size();
        for (uint i=0; i<n_indiv; i++)
        {
            if (_world[k].get_type() != SP_other)  // <-- do not change link while present in this SP
            {
                std::uniform_real_distribution<> unif01(0.0, 1.0);
                if(unif01(_RANDOM_GENERATOR) < p)
                {
                    // Old link to be removed:
                    ID id_sp_other_old = _world[k].get_indiv(i).get_id_sp_other();
                    ID id_indiv        = _world[k].get_indiv(i).get_id();
                    _world[id_sp_other_old].remove_linked_indiv(id_indiv);

                    // Note: nothing to remove for the individual
                    //  bc the 'set_sp_other_link' will overwrite
                    // individual's linked 'other' social place.
                    
                    // New link to be added:
                    socialPlace *sp_other_new = pick_rnd_sp_other();
                    set_sp_other_link(k, i, *sp_other_new);
                }
            }
        }
    }
}











