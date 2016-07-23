//
//  simulation.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-06.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "simulation.h"


void Simulation::base_constructor(){
    
    _n_E = 0;
    _n_Ia = 0;
    _n_Is = 0;
    _n_R  = 0;
    _n_H  = 0;
    _incidence = 0;
    _horizon = -999;
    
    _ts_times.clear();
    _ts_E.clear();
    _ts_Ia.clear();
    _ts_Is.clear();
    _ts_R.clear();
    _ts_census_by_SP.clear();
    
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
    
    // individuals
    uint num_indiv			= (uint)(n_indiv * 3); // <-- make sure it's large enough
    vector<individual> many_indiv	= build_individuals(num_indiv,
                                                        sched,
                                                        dol_distrib,
                                                        doi_distrib);
    
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
    
    // individuals
    uint num_indiv			= (uint)(n_indiv*2); // <-- make sure it's large enough
    vector<individual> many_indiv	= build_individuals(num_indiv,
                                                        sched,
                                                        dol_distrib,
                                                        doi_distrib);
    
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
    
    // Schedule
    vector<double> timeslice {8.0/24, 16.0/24}; // must sum up to 1.0
    vector<SPtype> worker_sed  {SP_workplace, SP_household};
    
    schedule sched_worker_sed(worker_sed, timeslice, "worker_sed");
    
    vector<schedule> sched {
        sched_worker_sed,
    };
    
    // Disease stage durations:
    string dol_distrib = "exp";
    string doi_distrib = "exp";
    
    // individuals
    uint num_indiv			= (uint)(n_indiv*2); // <-- make sure it's large enough
    vector<individual> many_indiv	= build_individuals(num_indiv,
                                                        sched,
                                                        dol_distrib,
                                                        doi_distrib);
    
    // type of social places
    // existing in the simulated world:
    vector<SPtype> spt {
        SP_workplace,
        SP_household,
        SP_hospital
    };
    
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
    
    // individuals
    uint num_indiv			= (uint)(1e7 * sizereduction); // <-- make sure it's large enough
    vector<individual> many_indiv	= build_individuals(num_indiv,
                                                        sched,
                                                        dol_distrib,
                                                        doi_distrib);
    
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
    
    
    // ----- MAIN LOOP FOR TIME ------
    
    if(debug_mode) check_book_keeping();
    
    define_all_id_tables();
    
    for (_current_time=0.0; _current_time < _horizon; ) {
        
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
        
        // Record for time series:
        _ts_times.push_back(_current_time);
        _ts_incidence.push_back(_incidence);
        _ts_prevalence.push_back(_prevalence);
        _ts_S.push_back(_n_S);
        _ts_E.push_back(_n_E);
        _ts_Ia.push_back(_n_Ia);
        _ts_Is.push_back(_n_Is);
        _ts_R.push_back(_n_R);
        
        // Advance time:
        time_update(dt);
        update_pop_count();
        k++;
        
        if(debug_mode) check_book_keeping();
    }
    
    if(debug_mode){
        cout << endl << endl << "Simulation completed."<< endl;
        displayVector(_ts_incidence);
        cout <<"Final size: " << sumElements(_ts_incidence)<<endl;
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
    // add individual at destination
    _world[to].add_indiv(tmp);
    // remove this individual (in pos_indiv^th position in '_indiv' vector) from here
    _world[from].remove_indiv(pos_indiv);
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
                    
                    if ( (id_dest!=k) && (id_dest!=__UNDEFINED_ID) )
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
    
    vector<uint> n_contacts; // number of contacts for each infectious individual
    
    uint n=0;
    if(infectious_type == "Is") n = (uint)_world[k]._indiv_Is.size();
    if(infectious_type == "Ia") n = (uint)_world[k]._indiv_Ia.size();
    
    for (uint i=0; i<n; i++) {
        // draw the contact rate
        // based on the individual and the current social place:
        individual* tmp = nullptr;
        if(infectious_type == "Is") tmp = _world[k]._indiv_Is[i];
        if(infectious_type == "Ia") tmp = _world[k]._indiv_Ia[i];
        
        double cr = draw_contact_rate(tmp, k);
        
        // Draw the actual number of contacts
        std::poisson_distribution<> poiss(cr * dt);
        n_contacts.push_back( poiss(_RANDOM_GENERATOR) );
    }
    return n_contacts;
}


double Simulation::calc_proba_transmission(individual *infectious,
                                           individual *susceptible){
    /// Calculate probability of transmission given contact
    /// between an infectious and susceptible individuals
    
    // TO DO: implement something more elaborate!
    
    // Susceptible side:
    double p_susc = susceptible->get_frailty() * (1 - susceptible->get_immunity());
    
    // infectious side:
    double p_inf = 1.0;
    if(!infectious->is_symptomatic()) p_inf = _modelParam.get_prm_double("asymptom_infectiousness_ratio");
    
    return p_susc * p_inf;
}


vector< vector<uint> > Simulation::draw_contacted_S(uint k,
                                                    vector<uint> n_contacts,
                                                    string infectious_type){
    /// Randomly select the susceptible that will be
    /// contacted by infectious individuals of a
    /// given infectious type (e.g. symptomatic or not).
    
    uint n_inf = 0;
    if(infectious_type == "Is") n_inf = (uint)_world[k]._indiv_Is.size();
    if(infectious_type == "Ia") n_inf = (uint)_world[k]._indiv_Ia.size();
    
    uint n_susc = (uint)_world[k]._indiv_S.size();
    vector< vector<uint> > selected_S(n_inf);
    
    // For all infectious indiv,
    // sample (with replacement)
    // which susceptibles will be contacted:
    std::uniform_int_distribution<uint> unif(0,n_susc-1);

    for(uint i=0; i<n_inf; i++){
        for(uint j=0; j<n_contacts[i]; j++)
        {
            uint rnd_idx = unif(_RANDOM_GENERATOR);
            selected_S[i].push_back(rnd_idx);
        }
    }
    return selected_S;
}


vector< vector<uint> > Simulation::transmission_attempts(uint k,
                                                         vector< vector<uint> > selected_S){
    /// Attempts transmission on all selected susceptible individuals
    
    uint n_inf = (uint) selected_S.size();
    vector< vector<uint> > transm_success(n_inf);
    
    bool homog = _modelParam.get_prm_bool("homogeneous_contact");
    
    if(homog){
        // Homogeneous contact: for testing purpose only. (but keep it!)
        // Always success because the contact rate
        // is understood as an _effective_ one in the homogeneous case:
         for (uint i=0; i<n_inf; i++)
             transm_success[i].resize(selected_S[i].size(), true);
    }
    
    if(!homog){
        std::uniform_real_distribution<float> unif01(0.0, 1.0);
        
        for (uint i=0; i<n_inf; i++)
        {
            transm_success[i].resize(selected_S[i].size(), false);
            
            for (uint j=0; j < selected_S[i].size(); j++){
                // TO DO: could be optimized:
                // because sampling with replacement,
                // we may attempt transmission on a susceptible
                // who had successful transmission previously
                
                double p_ij = calc_proba_transmission(_world[k]._indiv_Is[i],
                                                      _world[k]._indiv_S[selected_S[i][j]]);
                if (unif01(_RANDOM_GENERATOR)<p_ij) {
                    transm_success[i][j] = true;
                }
            }
        }
    }

    return transm_success;
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
    
    for (uint m=0; m<SS_melt_success.size(); m++) {
        _world[k]._indiv_S[SS_melt_success[m]]->acquireDisease();
        
        // Book keeping: Update counters
        _world[k].update_epidemic_count(*_world[k]._indiv_S[SS_melt_success[m]], "new_case");
        
        // Book keeping: Update pointer tables:
        uint tmp1 = SS_melt_success[m];
        _world[k]._indiv_S.erase(_world[k]._indiv_S.begin() + tmp1);
        
        incidence ++;
        // TO DO: record infection times, etc.
    }
    
    
//    for (uint i=0; i<selected_S.size(); i++) {
//        for (uint j=0; j < selected_S[i].size(); j++){
//            
//            if(transm_success[i][j]){
//                _world[k]._indiv_S[selected_S[i][j]]->acquireDisease();
//                
//                // Book keeping: Update counters
//                _world[k].update_epidemic_count(*_world[k]._indiv_S[selected_S[i][j]], "new_case");
//
//                // Book keeping: Update pointer tables:
//                uint tmp1 = selected_S[i][j];
//                _world[k]._indiv_S.erase(_world[k]._indiv_S.begin() + selected_S[i][j]);
//                
//                
//                incidence ++;
//                // TO DO: record infection times, etc.
//            }
//        }
//    }

    return incidence;
}


uint Simulation::transmission_oneSP(uint k,
											double contact_rate,
											double dt){
	/// Performs transmission within the k^th social place.
	/// Returns incidence for _this_ social place, during the time step 'dt'.
	
    // STOPPED HERE & TO DO:
    // implement the homogeneous contact case
    
	stopif(k >= _world.size(), "Asking for an inexistent social place");

    // If there are no susceptibles
    // or no infectious individuals,
    // then no transmission can occur!!!
    if (_n_S == 0 || (_n_Is+_n_Ia==0) ) {
        return 0;
    }
    
    // Number of contacts for each infectious individual.
    // (numbers are stored in the order of vectors _indiv_Is, _indiv_Ia)
    
    vector<uint> n_contacts_Is = draw_n_contacts(k, dt, "Is");
    vector<uint> n_contacts_Ia = draw_n_contacts(k, dt, "Ia");
    
    // Randomly select susceptibles
    // in this social place that will
    // be in contact: sampling _with_ replacement
    // (a S can be in contact with more than one I)
    
    vector< vector<uint> > selected_S_Is = draw_contacted_S(k, n_contacts_Is, "Is");
    vector< vector<uint> > selected_S_Ia = draw_contacted_S(k, n_contacts_Ia, "Ia");
    
    // Transmission attempts
    // (need to do this intermediary step
    // in order not to mess up pointers vector '_indiv_X')
    
    vector< vector<uint> > transm_success_Is = transmission_attempts(k, selected_S_Is);
    vector< vector<uint> > transm_success_Ia = transmission_attempts(k, selected_S_Is);
    
    // Activate disease acquisition:
    
    uint inc_Is     = transmission_activation(k, selected_S_Is, transm_success_Is);
    uint inc_Ia     = transmission_activation(k, selected_S_Ia, transm_success_Ia);
    uint inc_total  = inc_Is + inc_Ia;
    
	return inc_total;
}


void Simulation::transmission_world(double timeslice){
    /// Simulates disease transmissions in the whole world (all social places)
    
    uint incidence = 0;
    double cr = _modelParam.get_prm_double("contact_rate");
    
    for(uint k=0; k < _world.size(); k++){
        incidence += transmission_oneSP(k, cr, timeslice);
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
    
    ID cnt = 0;
    for(ID i=0; i<_world.size(); i++){
        if (_world[i].get_id_sp() == id_sp[cnt]) {
            
            // check:
            string errmsg =  "Cannot seed in SP_ID_" + to_string(i) + " because it has less individuals than intial infections requested!";
            stopif(_world[i].get_size() < I0[cnt], errmsg);
            
            // seed infection in this social place:
            for(uint k=0; k<I0[cnt]; k++)	_world[i].acquireDisease(k);
            _world[i].set_n_E(I0[cnt]);
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
    
}



void Simulation::define_all_id_tables(){
    /// Define all IDs of susceptible
    /// and infectious individuals
    /// in all social places
    /// (define _id_S for all social places)
    
    for (uint k=0; k<_world.size(); k++)
    {
        vector<individual> indiv = _world[k].get_indiv();

        _world[k].clear_id_S();
        _world[k].clear_id_Is();
        _world[k].clear_id_Ia();
        _world[k]._indiv_S.clear();
        _world[k]._indiv_Is.clear();
        _world[k]._indiv_Ia.clear();
        
        for (uint i=0; i< indiv.size(); i++)
        {
            if (indiv[i].is_susceptible()) {
                _world[k].add_id_S(indiv[i].get_id());
                _world[k]._indiv_S.push_back(_world[k].get_mem_indiv(i));
            }
            else if (indiv[i].is_infectious() && indiv[i].is_symptomatic()){
                _world[k].add_id_Is(indiv[i].get_id());
                _world[k]._indiv_Is.push_back(_world[k].get_mem_indiv(i));
            }
            else if (indiv[i].is_infectious() && !indiv[i].is_symptomatic()){
                _world[k].add_id_Ia(indiv[i].get_id());
                _world[k]._indiv_Ia.push_back(_world[k].get_mem_indiv(i));
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


// STOPPED HERE: change the way hospitalization is done
// Instead of checking at every time step among Is,
// decide when infection is acquired:
// - if indiv will be hospitalized
// - when
// - for how long (doh)
// But then, must implement a check if for example
// treatment is given before hospitalization,
// hospitalization must be canceled and indiv moved to T compartment.


void Simulation::hospitalize(){
    /// Scan all infectious symptomatic individuals
    /// and 'decide' (stocahstic decision) to hospitalize or not
    
    std::uniform_real_distribution<> unif(0.0, 1.0);
    
    for (uint k=0; k<_world.size(); k++) {
        uint nIs = _world[k].get_n_Is();
        for (uint i=0; i<nIs; i++) {
            
            // If not already hospitalized:
            if (! _world[k]._indiv_Is[i]->is_hosp())
            {
                double p =_modelParam.get_prm_double("proba_hospitalization");
                if (unif(_RANDOM_GENERATOR) < p)
                {
                    // set the hospitalization flag for this individual
                    _world[k]._indiv_Is[i]->set_is_hosp(true);
                    
                    // Move individual to its linked SP_hospital
                    ID id_sp_hospital = _world[k]._indiv_Is[i]->get_id_sp_hospital();
                    ID id_indiv = _world[k]._indiv_Is[i]->get_id();
                    uint pos = _world[k].find_indiv_pos(id_indiv);
                    move_one_individual(pos, k, id_sp_hospital);
                    
                    // update counters
                    _n_H++;
                    _world[id_sp_hospital].increment_n_H();
                    
                    // update _indiv_Is pointers:
                    _world[id_sp_hospital]._indiv_Is.push_back(_world[k]._indiv_Is[i]);  //_world[k].get_mem_indiv(pos)
                    _world[k]._indiv_Is.erase(_world[k]._indiv_Is.begin()+i);
                }
            }
        }
    }
}













