//
//  simulation.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-06.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "simulation.h"


void Simulation::base_constructor(){
    
    _n_E  = 0;
    _n_Ia = 0;
    _n_Is = 0;
    _n_R  = 0;
    _n_H  = 0;
    _n_D  = 0;
    _incidence = 0;
    _n_treated = 0;
    _n_vaccinated = 0;
    _horizon = -999;
    
    _ts_times.clear();
    _ts_E.clear();
    _ts_Ia.clear();
    _ts_Is.clear();
    _ts_R.clear();
    _ts_D.clear();
    _ts_census_by_SP.clear();
    _ts_n_treated.clear();
    _ts_n_vaccinated.clear();
    
    _intervention.clear();
}


Simulation::Simulation(){
    base_constructor();
}


void Simulation::build_single_world(uint n_indiv){
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
    probaDistrib<uint> p_hh({n_indiv},{1.0});
    
    vector< probaDistrib<uint> > p_size {p_hh};
    
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


void Simulation::build_test_2_sp(uint n_indiv){
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
    probaDistrib<uint> p_workPlace({n_indiv/2},  {1.0});
    probaDistrib<uint> p_hh({n_indiv},  {1.0});
    
    vector<probaDistrib<uint> > p_size {
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



void Simulation::build_test_hospitalization(uint n_indiv){
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
    probaDistrib<uint> p_workPlace({n_indiv/2},  {1.0});
    probaDistrib<uint> p_hh({n_indiv},  {1.0});
    probaDistrib<uint> p_hosp({0},  {1.0});
    
    vector<probaDistrib<uint> > p_size {
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




void Simulation::build_test_world(double sizereduction){
    
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
    
    probaDistrib<uint> p_pubTransp({20,30,60},{0.6,0.3,0.1});
    probaDistrib<uint> p_workPlace({3,7,15,30,75,200},{0.6,0.15,0.15,0.07,0.02,0.01});
    probaDistrib<uint> p_hh({1,2,3,4,5,6,7,8},{0.23, 0.34, 0.16, 0.15, 0.06, 0.03, 0.02, 0.01});
    probaDistrib<uint> p_school({250,500,750,1000,1250,1500},{0.60,0.25,0.10,0.03, 0.01,0.01});
    
    vector<probaDistrib<uint> > p_size {
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


void Simulation::time_update(double dt){
    /// Make all relevant updates when time is advanced by 'dt'
    
    // Main timer:
    _current_time += dt;
    
    // Update individuals' clock:
    for (uint k=0; k<_world.size(); k++) {
        _world[k].time_update(dt);
    }
    update_pop_count();
}


void Simulation::test(){
    
    //    _world[0]._indiv_S[0]->set_immunity(0.12345);
    
}


void Simulation::run(){
    /// Run the simulated epidemic
    
    
    // Retrieve all model parameters:
    double p_move = _modelParam.get_prm_double("proba_move");
    bool debug_mode = _modelParam.get_prm_bool("debug_mode");
    
    if(debug_mode) cout << endl << endl << " ======= START SIMULATION ======" <<endl<<endl;
    
    _current_time = 0.0;
    
    // TO DO: CHANGE THAT, IT's UGLY AND DANGEROUS
    ID ii = at_least_one_indiv_present(_world)[0];
    vector<double> timeslice = _world[ii].get_indiv()[0].get_schedule().get_timeslice();
    // - - - - - - - - -
    
    unsigned long nts = timeslice.size();
    uint k = 0;
    
    if(debug_mode) check_book_keeping();
    define_all_id_tables();
    
    
    // ----- MAIN LOOP FOR TIME ------
    
    for (_current_time = 0.0; _current_time < _horizon; ) {
        
        uint idx_timeslice = k % nts;
        double dt = timeslice[idx_timeslice];
        
        if(debug_mode){
            cout << "iter = " << k;
            cout << " ; time = " << _current_time;
            cout << " ; sched idx = " << idx_timeslice;
            cout << " ; sched length = " << dt;
            cout <<endl;
        }
        
        update_ts_census_by_SP();
        
        discharge_hospital(idx_timeslice);
        define_all_id_tables(); // <-- check if this is necessary here
        
        death_hospital();
        define_all_id_tables(); // <-- check if this is necessary here
        
        if(p_move>0) {
            move_individuals_sched(idx_timeslice, p_move);
            define_all_id_tables();
            if(debug_mode) check_book_keeping();
        }
        
        transmission_world(dt);
        update_pop_count();
        
        if(debug_mode) check_book_keeping();
        
        hospitalize();
        define_all_id_tables(); // <-- check if this is necessary here
        update_pop_count(); // <-- check if this is necessary here
        
        if(debug_mode) check_book_keeping();
        
        // Interventions
        for (uint k=0; k<_world.size(); k++) {
            activate_interventions(k, dt);
        }
        
        if(debug_mode) check_book_keeping();
        
        if(at_least_one_vaccination_intervention())
            update_immunity_frailty();
        
        // Record for time series:
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
        _ts_n_treated.push_back(_n_treated);
        _ts_n_vaccinated.push_back(_n_vaccinated);
        
        // Advance time:
        time_update(dt);
        k++;
        
        if(debug_mode) check_book_keeping();
    }
    
    if(debug_mode){
        cout << endl << endl << "Simulation completed."<< endl;
        
        cout << "time series of incidence:"<<endl;
        displayVector(_ts_incidence);
        
        cout <<"Final size: " << sumElements(_ts_incidence)<<endl;
        
        cout << endl << "time series of hospitalizations:"<<endl;
        displayVector(_ts_H);
        
        cout << endl << "time series of deaths:"<<endl;
        displayVector(_ts_D);
    }
}


void Simulation::set_world(world w){
    _world = w;
    update_pop_count();
}

void Simulation::set_disease(const disease &d){
    /// Set the disease 'd' to all individuals in all social places
    for (uint k=0; k<_world.size(); k++) {
        _world[k].set_disease_to_all_indiv(d);
    }
}


void Simulation::move_one_individual(uint pos_indiv, ID from, ID to){
    /// Move the individual in position "pos_indiv" in thevector "_indiv"
    /// from one social place to another.
    /// (social places are identified by their IDs/position)
    
    individual tmp = _world[from].get_indiv(pos_indiv);
    
    // DEBUG
    ID tmpid =_world[from].get_indiv(pos_indiv).get_id();
    
    
    // add individual at destination
    _world[to].add_indiv(tmp);
    // remove this individual (in pos_indiv^th position in '_indiv' vector) from here
    _world[from].remove_indiv(pos_indiv);
    
    // DEBUG
    //cout << "DEBUG --> id_" <<tmpid << "moved from SP_"<<from << " to SP_"<<to <<endl;
}


void Simulation::move_individuals_sched(uint idx_timeslice,
                                        double proba){
    /// Move individuals across social places according to their schedule
    
    unsigned long N = _world.size();
    
    std::uniform_real_distribution<double> unif(0.0,1.0);
    
    for (int k=0; k<N; k++)
    {
        uint n = (uint)_world[k].get_size();
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
                    
                    stopif(id_dest==__UNDEFINED_ID,
                           "Undefined destination for indiv ID_" + to_string(_world[k].get_indiv(i).get_id()));
                    
                    if ( (id_dest!=k) )
                    {
                        // take a copy of the individual
                        individual tmp = _world[k].get_indiv(i);
                        
                        // Draw the chance move will actually happen:
                        // TO DO: make this proba individual-dependent
                        double u = unif(_RANDOM_GENERATOR);
                        if ( u < proba ){
                            
                            move_one_individual(i, k, id_dest);
                            
                            // DELETE WHEN SURE:
                            //                        // add individual at destination
                            //                        _world[id_dest].add_indiv(tmp);
                            //                        // remove this individual (in i^th position in '_indiv' vector) from here
                            //                        _world[k].remove_indiv(i);
                        }
                    }
                } // end-if-not-hospitalized
            }
        }
    } // end-for-k-socialPlace
} // end-function


void Simulation::move_individuals(const SPtype sptype, double proba){
    
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


vector<uint> Simulation::draw_n_contacts(uint k,
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


double Simulation::calc_proba_transmission(individual *infectious,
                                           individual *susceptible){
    /// Calculate probability of transmission given contact
    /// between an infectious and susceptible individuals
    
    // TO DO: implement something more elaborate!
    
    // Susceptible side (acquisition risk):
    double p_susc = susceptible->get_frailty() * (1 - susceptible->get_immunity());
    
    // infectious side (transmission risk):
    double p_inf = 1.0;
    if(! infectious->is_symptomatic())
        p_inf = _modelParam.get_prm_double("asymptom_infectiousness_ratio");
    
    if (infectious->is_treated()){
        // Infectiousness reduction if patient is treated
        //
        double m = _modelParam.get_prm_double("treat_reduc_infect_mean");
        if(m>0){
            double reduc = beta_distribution(1.0, 1/m - 1.0, _RANDOM_GENERATOR);
            p_inf = p_inf * (1-reduc);
        }
    }
    
    // DEBUG
//    cout << " DEBUG: p_susc = "<< p_susc << " ; p_inf = " << p_inf << endl;
    
    return p_susc * p_inf;
}



double Simulation::calc_proba_symptomatic(float immunity, float frailty){
    /// Probability to be symptomatic given
    /// an individual's immunity and frailty
    
    // TO DO: more sophisticated!
    double res = frailty * (1-immunity); //(frailty * (1-immunity) < 0.5)? 0.2 : 1.0;
    return res;
}


double Simulation::calc_proba_hospitalized(float frailty){
    /// Probability to be symptomatic given
    /// an individual's frailty
    
    // TO DO: more sophisticated!
    
    return frailty ;
}


double Simulation::calc_proba_death(float frailty){
    /// Probability to die at the end of hospitalization period.
    
    // TO DO: more sophisticated!
    double p        = _modelParam.get_prm_double("proba_death_prm_1");
    double thres    = _modelParam.get_prm_double("proba_death_prm_2");
    double p_hi     = _modelParam.get_prm_double("proba_death_prm_3");
    if (frailty > thres) p = p_hi;
    return p;
}



vector< vector<uint> > Simulation::draw_contacted_S(uint k,
                                                    vector<uint> n_contacts,
                                                    string infectious_type){
    /// Randomly select the positions in '_indiv_S'
    /// for susceptibles that will be
    /// contacted by infectious individuals of a
    /// given infectious type (e.g. symptomatic or not).
    
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


vector< vector<uint> > Simulation::transmission_attempts(uint k,
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


void Simulation::transmission_wiw(int k,
                                  vector<vector<uint> > selected_S,
                                  vector<vector<uint> > transm_success,
                                  string infectious_type){
    /// Records who infected who.
    
    for (uint i=0; i<selected_S.size(); i++) {
        for (uint s=0; s<selected_S[i].size(); s++) {
            // if transmission is successful:
            if (transm_success[i][s])
            {
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
            } // end-if-transmission-success
        } //end-for-s
    } //end-for-i
}

uint Simulation::transmission_activation(int k,
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
        float immunity = _world[k]._indiv_S[SS_melt_success[m]]->get_immunity();
        float frailty  = _world[k]._indiv_S[SS_melt_success[m]]->get_frailty();
        double p_sympt = calc_proba_symptomatic(immunity, frailty);
        
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



uint Simulation::transmission_process(uint k, double dt, string infectious_type){
    /// Full transmission process:
    /// draw number of contacts, identify susceptible contacted,
    /// attempts transmission, activate successful attempts.
    
    stopif(k >= _world.size(), "Asking for an inexistent social place");
    
    // If there are no susceptibles
    // or no infectious individuals,
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
    vector< vector<uint> > selected_S_Ix = draw_contacted_S(k, n_contacts_Ix, infectious_type);

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


uint Simulation::transmission_oneSP(uint k,
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


void Simulation::transmission_world(double timeslice){
    /// Simulates disease transmissions in the whole world (all social places)
    
    uint incidence = 0;
    for(uint k=0; k < _world.size(); k++){
        incidence += transmission_oneSP(k, timeslice);
    }
    _incidence = incidence;
}


uint Simulation::census_total_alive(){
    /// Counts all individuals that are alive
    uint cnt = 0;
    for(int k=0; k<_world.size(); k++) cnt += _world[k].census_alive();
    return cnt;
}


uint Simulation::prevalence(){
    /// Prevalence in the whole world
    
    uint cnt = 0;
    for (int k=0; k<_world.size(); k++) cnt += _world[k].get_prevalence();
    return cnt;
}


uint Simulation::population_size(){
    
    uint s = 0;
    for(int i=0; i<_world.size(); i++) s+=_world[i].get_size();
    return s;
}


void Simulation::display_split_pop_present(){
    
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

void Simulation::display_split_pop_linked(){
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

void Simulation::seed_infection(vector<ID> id_sp, vector<uint> I0){
    /// Seed infection in specified socialplaces, with specified initial number of infectious indiv
    
    stopif(id_sp.size() != I0.size(), "vectors must be same size");
    
    vector<vector<uint> > selected_S(1);
    vector<vector<uint> > transm_success(id_sp.size());
    
    ID cnt = 0;
    for(ID i=0; i<_world.size(); i++){
        if (_world[i].get_id_sp() == id_sp[cnt]) {
            
            // check:
            string errmsg =  "Cannot seed in SP_ID_" + to_string(i) + " because it has less individuals than intial infections requested!";
            stopif(_world[i].get_size() < I0[cnt], errmsg);
            
            // seed infection in this social place:
            //            for(uint k=0; k<I0[cnt]; k++)	_world[i].acquireDisease(k);
            //            _world[i].set_n_E(I0[cnt]);
            
            for(uint j=0; j<I0[cnt]; j++) {
                selected_S[0].push_back(j);
                transm_success[0].push_back(true);
            }
            transmission_activation(id_sp[cnt], selected_S, transm_success);
            
            cnt++;
        }
    }
    update_pop_count();
}


void Simulation::displayInfo_indiv(){
    /// Display informations on all individuals, in all social places:
    
    for (uint k=0; k<_world.size(); k++) {
        for (ID i=0; i<_world[k].get_size(); i++) {
            _world[k].get_indiv(k).displayInfo();
        }
    }
}


dcDataFrame Simulation::timeseries(){
    
    dcDataFrame df(_ts_times,"time");
    
    df.addcol("incidence", to_vector_double(_ts_incidence));
    df.addcol("prevalence", to_vector_double(_ts_prevalence));
    df.addcol("nS", to_vector_double(_ts_S));
    df.addcol("nE", to_vector_double(_ts_E));
    df.addcol("nIa", to_vector_double(_ts_Ia));
    df.addcol("nIs", to_vector_double(_ts_Is));
    df.addcol("nR", to_vector_double(_ts_R));
    df.addcol("n_treated", to_vector_double(_ts_n_treated));
    df.addcol("n_vaccinated", to_vector_double(_ts_n_vaccinated));
    
    return df;
}



void Simulation::update_ts_census_by_SP(){
    /// Update the data frame recording
    /// the number of individuals in each
    /// disease stage, for every social places.
    
    unsigned long n = _world.size();
    
    if (_ts_census_by_SP.get_nrows() == 0)
    {
        // Create data frame with time
        vector<double> t(n, 0.0);
        dcDataFrame tmp(t,"time");
        _ts_census_by_SP = tmp;
        
        vector<double> id_sp(n);
        vector<double> pop_present(n);
        vector<double> nS(n);
        vector<double> nE(n);
        vector<double> nIs(n);
        vector<double> nIa(n);
        vector<double> nR(n);
        
        for (uint k=0; k<n; k++) {
            id_sp[k]        = _world[k].get_id_sp();
            pop_present[k]  = _world[k].get_size();
            nS[k]           = _world[k].get_n_S();
            nE[k]           = _world[k].get_n_E();
            nIs[k]          = _world[k].get_n_Is();
            nIa[k]          = _world[k].get_n_Ia();
            nR[k]           = _world[k].get_n_R();
        }
        
        _ts_census_by_SP.addcol("id_sp", id_sp);
        _ts_census_by_SP.addcol("pop_present", pop_present);
        _ts_census_by_SP.addcol("nS", nS);
        _ts_census_by_SP.addcol("nE", nE);
        _ts_census_by_SP.addcol("nIs", nIs);
        _ts_census_by_SP.addcol("nIa", nIa);
        _ts_census_by_SP.addcol("nR", nR);
    }
    
    if (_ts_census_by_SP.get_nrows() > 0)
    {
        // Create data frame with time
        vector<double> t(n, _current_time);
        dcDataFrame tmp(t,"time");
        
        vector<double> id_sp(n);
        vector<double> pop_present(n);
        vector<double> nS(n);
        vector<double> nE(n);
        vector<double> nIs(n);
        vector<double> nIa(n);
        vector<double> nR(n);
        
        for (uint k=0; k<n; k++) {
            id_sp[k]        = _world[k].get_id_sp();
            pop_present[k]  = _world[k].get_size();
            nS[k]           = _world[k].get_n_S();
            nE[k]           = _world[k].get_n_E();
            nIs[k]          = _world[k].get_n_Is();
            nIa[k]          = _world[k].get_n_Ia();
            nR[k]           = _world[k].get_n_R();
        }
        
        tmp.addcol("id_sp", id_sp);
        tmp.addcol("pop_present", pop_present);
        tmp.addcol("nS", nS);
        tmp.addcol("nE", nE);
        tmp.addcol("nIs", nIs);
        tmp.addcol("nIa", nIa);
        tmp.addcol("nR", nR);
        
        _ts_census_by_SP = rbind(_ts_census_by_SP, tmp);
    }
    
    // DEBUG
    //    _ts_census_by_SP.display();
}


void Simulation::update_pop_count(){
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


void Simulation::check_book_keeping(){
    /// Check consistency of book keeping.
    /// WARNING: FOR DEBUG ONLY, SLOWS DOWN EXECUTION!
    
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



void Simulation::define_all_id_tables(){
    /// Define all IDs and pointers
    /// of tracked individuals
    /// in all social places.
    
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


double  Simulation::draw_contact_rate(individual* indiv, uint k){
    /// Draw the contact rate for an infectious individual
    /// in a given social place.
    
    bool homog = _modelParam.get_prm_bool("homogeneous_contact");
    
    double cr = -999.99;
    
    if(homog){
        // Homogeneous contacts.
        // Used for testing purposes (e.g. compare to ODE models).
        uint ns = _world[k].get_n_S();
        uint n  = (uint)_world[k].get_size();
        cr = _modelParam.get_prm_double("contact_rate") * ns / n;
    }
    
    if(!homog){
        // TO DO: implement something better. Test for now...
        double mult_indiv = 1.0;
        double mult_sp = 1.0;
        SPtype sp_type = _world[k].get_type();
        
        double age = indiv->get_age();
        if (age <30.0) mult_indiv = 2.0;
        if (sp_type == SP_household) mult_sp = 1.9;
        if (sp_type == SP_pubTransp) mult_sp = 1.5;
        
        cr = mult_indiv * mult_sp * _modelParam.get_prm_double("contact_rate");
    }
    return cr;
}


//void Simulation::hospitalize_indiv(uint k, uint i){
//    /// Hospitalize individual #i on social place #k
//
//}

void Simulation::assign_hospital_to_individuals(){
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


void Simulation::hospitalize(){
    /// Scan all infectious symptomatic individuals
    /// and hospitalize them if they were meant to and if it's time to do so.
    
    std::uniform_real_distribution<> unif(0.0, 1.0);
    
    for (uint k=0; k<_world.size(); k++)
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
                
                // update counters
                _n_H++;
                _world[id_sp_hospital].increment_n_H();
                
                
                // Move individual to its linked SP_hospital
                move_one_individual(pos_indiv, k, id_sp_hospital);
//                cout <<" DEBUG: indiv ID_"<<to_string(id_indiv)<<"  hospitalized" <<endl;
            }
        }
    }
}


void Simulation::discharge_hospital(uint idx_timeslice){
    /// Scan all infectious hospitalized individuals
    /// and discharge them if duration of hospitalization > than the one drawn (_doh_drawn)
    
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
                    
                    // DELETE WHEN SURE:
//                    uint pos_Is = _world[k].find_indiv_X_pos(id_indiv, "Is");
//                    _world[k]._indiv_Is[pos_Is]->set_is_hosp(false);

                    
                    // update counters
                    _n_H--;
                    _world[id_sp_hospital].decrement_n_H();
                    
                    // Move individual to its scheduled social place.
                    // Retrieve its actual destination:
                    ID id_dest = _world[k].find_dest(i, idx_timeslice);
                    
                    stopif(id_dest==__UNDEFINED_ID,
                           "Undefined destination for indiv ID_" + to_string(id_indiv));
                    
                    if ( id_dest!=k )
                    { move_one_individual(pos_indiv, k, id_dest); }
                    
                    //                    cout <<" DEBUG: indiv ID_"<<to_string(id_indiv)<<" leave hospital!" <<endl;
                }
            }
        } // end-if-SP_hospital
    }
}

void Simulation::death_hospital(){
    
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




vector<individual*> Simulation::draw_targeted_individuals(uint i,
                                                           ID id_sp,
                                                           double dt){
    /// Draw the targeted individuals of the ith intervention,
    /// in sociale place with ID 'id_sp'
    
    vector<individual*> indiv_drawn;
    bool found = false;
    
    float cvg_rate  = _intervention[i].get_cvg_rate();
    string type_target = _intervention[i].get_type_indiv_targeted();
    float interv_intensity = cvg_rate * _world[id_sp].get_size() * dt;
    
    if(type_target == "symptomatic"){
        
        // * * WARNING * *  excludes symptomatic already treated!
        
        found = true;
        uint nI = (uint)_world[id_sp]._indiv_Is.size();
        
        if (nI > 0){
            // Draw the total number of individuals that are targeted
            std::poisson_distribution<> poiss(interv_intensity);
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
        
        // * * WARNING * * excludes susceptiblr already treated & vaccinated!
        
        found = true;
        uint nS  = (uint)_world[id_sp]._indiv_S.size();

        if (nS > 0){
            // Draw the total number of individuals that are targeted
            std::poisson_distribution<> poiss(interv_intensity);
            uint n_target = poiss(_RANDOM_GENERATOR);
            if (n_target > nS) n_target = nS;
            // Select targeted individuals:
            uint cnt = 0;
            for(uint i=0; i<nS && cnt<n_target; i++) {
                bool is_treated = _world[id_sp]._indiv_S[i]->is_treated();
                bool is_vax     = _world[id_sp]._indiv_S[i]->is_vaccinated();
                if (!is_treated && !is_vax){
                    indiv_drawn.push_back(_world[id_sp]._indiv_S[i]);
                    cnt ++;
                }
            }
        }
    } // end-if-type_target == "susceptible"
    
    
    stopif(!found, "Type of targeted individual unknown: " + type_target);
    
    return indiv_drawn;
}


void Simulation::activate_interventions(ID id_sp, double dt){
    /// Activate all interventions for social place 'id_sp'
    
    float treat_doi_reduc   = _modelParam.get_prm_double("treat_doi_reduc");
    float vax_imm_incr      = _modelParam.get_prm_double("vax_imm_incr");
    float vax_frail_incr    = _modelParam.get_prm_double("vax_frail_incr");
    float vax_lag           = _modelParam.get_prm_double("vax_lag_full_efficacy");
    
    for (uint i=0; i<_intervention.size(); i++) {
        
        if (_intervention[i].get_time_start() <= _current_time &&
            _intervention[i].get_time_end() > _current_time)
        {
            vector<individual*> x = draw_targeted_individuals(i, id_sp, dt);
            _intervention[i].act_on_individual(x,
                                               _current_time,
                                               treat_doi_reduc,
                                               vax_imm_incr,
                                               vax_frail_incr,
                                               vax_lag);
            
            if (_intervention[i].get_type_intervention()=="treatment")
                _n_treated += x.size();
            else if (_intervention[i].get_type_intervention()=="vaccination"){
                _n_vaccinated += x.size();
                
                _world[id_sp]._indiv_vax.insert(_world[id_sp]._indiv_vax.end(),
                                                x.begin(),
                                                x.end());
            }
        }
    }
}

void Simulation::update_immunity_frailty() {
    /// Update immunity and frailty of vaccinated individuals
    /// at each time step. Because vaccine takes some time
    /// to reach its full efficacy, both immunity and frailty
    /// need to be updated during this period of time.
    
    for (uint k=0; k<_world.size(); k++) {
        unsigned long n_vax = _world[k]._indiv_vax.size();
        for (uint i=0; i<n_vax; i++)
        {
            float t_vax   = _world[k]._indiv_vax[i]->get_vax_time_received();
            float lag_vax = _world[k]._indiv_vax[i]->get_vax_lag_full_efficacy();
            
            if (_current_time < t_vax+lag_vax){
                
                float target_imm =  _world[k]._indiv_vax[i]->get_vax_target_immunity();
                float target_fra =  _world[k]._indiv_vax[i]->get_vax_target_frailty();
                float imm0       =  _world[k]._indiv_vax[i]->get_imm_when_recv_vax();
                float frail0     =  _world[k]._indiv_vax[i]->get_frail_when_recv_vax();
                
                // Exponential growth to target level:
                float tmp_i_b = log(target_imm/imm0) / lag_vax;
                float tmp_i_a = imm0 * exp(-tmp_i_b * t_vax);
                float tmp_i   = tmp_i_a * exp( tmp_i_b * _current_time);
                
                float tmp_f_b = log(target_fra/frail0) / lag_vax;
                float tmp_f_a = frail0 * exp(-tmp_f_b * t_vax);
                float tmp_f   = tmp_f_a * exp( tmp_f_b * _current_time);
                
                // update levels:
                _world[k]._indiv_vax[i]->set_immunity(tmp_i);
                _world[k]._indiv_vax[i]->set_frailty(tmp_f);
            }
        }
    }
}



bool Simulation::at_least_one_vaccination_intervention(){
    /// Check if there is at least one vaccination among all
    /// interventions defined for this Simulator
    
    for (uint i=0; i<_intervention.size(); i++) {
        if (_intervention[i].get_type_intervention() == "vaccination")
            return true;
    }
    return false;
}









