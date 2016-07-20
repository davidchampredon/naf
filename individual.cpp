//
//  individual.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-04.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "individual.h"
#include "socialPlace.h"

void individual::base_constructor(){
    /// Base constructor for individuals
    /// (will be called by all genuine constructors)
    _is_alive	= true;
    _id			= __UNDEFINED_ID;
    _age		= 0.0;
    
    _immunity	= 0.0;
    _frailty	= 1.0;
    
    // Social place:
    _id_sp_current		= __UNDEFINED_ID;
    _id_sp_household	= __UNDEFINED_ID;
    _id_sp_workplace	= __UNDEFINED_ID;
    _id_sp_school		= __UNDEFINED_ID;
    _id_sp_other		= __UNDEFINED_ID;
    _id_sp_hospital		= __UNDEFINED_ID;
    _id_sp_pubTransp	= __UNDEFINED_ID;
    
    // Clinical:
    _is_at_risk		= false;
    
    _is_susceptible = true;
    _is_infected	= false;
    _is_latent		= false;
    _is_infectious	= false;
    _is_recovered	= false;
    _is_symptomatic = false;
    _was_symptomatic = false;
    _is_hosp		= false;
    
    _doi	= 0.0;
    _dol	= 0.0;
    _doh	= 0.0;
    _doi_drawn	= 0.0;
    _dol_drawn	= 0.0;
    _doh_drawn	= 0.0;
    
    _dol_distrib = "dirac";
    _doi_distrib = "dirac";
}


individual::individual(){
    base_constructor();
}


individual::individual(ID id, float age){
    base_constructor();
    _id = id;
    _age = age;
}



string individual::_disease_status_update(double dt){
    /// Updates disease-related member variables
    /// for this individual when time advances by 'dt'
    /// Returns a string that indicates what
    /// disease stage change (if any) occured.
    /// This string will be used at the social place level
    /// to update the counts of individuals.
    
    string res = "NO_CHANGE";
    
    // DEBUG
    //	if(_id == 0){
    //		cout << "debug"<<endl;
    //	}
    // ------
    
    // Disease-related durations:
    if (_dol>0 && _dol <= _dol_drawn) _dol += dt;
    if (_doi>0 && _doi <= _doi_drawn) _doi += dt;
    if (_doh>0 && _doh <= _doh_drawn) _doh += dt;
    
    // Update disease stages
    
    // From 'E' to 'I'
    if (_dol > _dol_drawn && !_is_infectious && !_is_recovered ) {
        _is_latent		= false;
        _is_infectious	= true;
        _doi = SUPERTINY;
        res = "E_to_I";
    }
    // From 'I' to 'R'
    if (_doi > _doi_drawn && !_is_recovered ) {
        _is_infectious	= false;
        _is_infected	= false;
        _is_symptomatic	= false;
        _is_recovered	= true;
        res = "I_to_R";
    }
    
    if (_doh > _doh_drawn) {
        _is_hosp = false;
        res = "LEAVE_HOSPITAL";
    }
    
    return res;
}


string individual::time_update(double dt){
    /// Updates all relevant member variables
    /// for this individual when time advances by 'dt'
    
    _age += dt/365.0;
    
    return _disease_status_update(dt);
}


void individual::set_id_sp_household(socialPlace& sp){
    stopif(sp.get_type() != SP_household, "social space must be a household!");
    _id_sp_household = sp.get_id_sp();
    sp.add_linked_indiv(_id);
}

void individual::set_id_sp_workplace(socialPlace& sp){
    stopif(sp.get_type() != SP_workplace, "social space must be a workplace!");
    _id_sp_workplace = sp.get_id_sp();
    sp.add_linked_indiv(_id);
}

void individual::set_id_sp_school(socialPlace& sp){
    stopif(sp.get_type() != SP_school, "social space must be a school!");
    _id_sp_school = sp.get_id_sp();
    sp.add_linked_indiv(_id);
}

void individual::set_id_sp_other(socialPlace& sp){
    stopif(sp.get_type() != SP_other, "social space must be a other public space!");
    _id_sp_other = sp.get_id_sp();
    sp.add_linked_indiv(_id);
}

void individual::set_id_sp_hospital(socialPlace& sp){
    stopif(sp.get_type() != SP_hospital, "social space must be a hospital!");
    _id_sp_hospital = sp.get_id_sp();
    sp.add_linked_indiv(_id);
}

void individual::set_id_sp_pubTransp(socialPlace& sp){
    stopif(sp.get_type() != SP_pubTransp, "social space must be a public transportation!");
    _id_sp_pubTransp = sp.get_id_sp();
    sp.add_linked_indiv(_id);
}

void individual::set_id_sp(SPtype type, socialPlace& sp){
    if (type == SP_hospital)  set_id_sp_hospital(sp);
    if (type == SP_household) set_id_sp_household(sp);
    if (type == SP_school)    set_id_sp_school(sp);
    if (type == SP_pubTransp) set_id_sp_pubTransp(sp);
    if (type == SP_other)     set_id_sp_other(sp);
    if (type == SP_workplace) set_id_sp_workplace(sp);
}


void individual::displayInfo(){
    
    cout << endl << "-- " << endl;
    cout << "individual ID: " << _id << endl;
    cout << "individual alive: " << _is_alive << endl;
    cout << "individual age: " << _age << endl;
    cout << "individual infected: " << _is_infected << endl;
    cout << "individual DOI: " << _doi << endl;
    cout << "individual's current SP id: " << _id_sp_current << endl;
    cout << "individual's household SP id: " << _id_sp_household << endl;
    cout << "individual's workplace SP id: " << _id_sp_workplace << endl;
    cout << "individual's school SP id: " << _id_sp_school << endl;
    cout << "individual's other pub. space SP id: " << _id_sp_other << endl;
    cout << "individual's hospital SP id: " << _id_sp_hospital << endl;
    cout << "individual's pubTransp SP id: " << _id_sp_pubTransp << endl;
    cout << "individual's schedule name: " << _schedule.get_name() << endl;
    cout << "name of the defined disease: " << _disease.get_name() << endl;
    cout << "dol = " << _dol << endl;
    cout << "doi = " << _doi << endl;
    cout << "doh = " << _doh << endl;
    cout << "-- " << endl;
}


double individual::calc_proba_acquire_disease(){
    /// Calculate the probability of aquiring
    /// the disease, given an infectious contact
    
    // TO DO: maybe implement something more sophisticated ? (logistic curves?)
    
    return (1-_immunity) * _frailty ;
}

void individual::acquireDisease(){
    /// This individual acquires the disease
    
    _is_susceptible = false;
    _is_infected	= true;
    _is_latent		= true;
    _is_infectious	= false;
    
    _dol = SUPERTINY;
    
    // Draws the durations THIS individual
    // will experience:
    
    // TO DO: implement something more sophisticated!
    
    if (_dol_distrib == "dirac") _dol_drawn = _disease.get_dol_mean();
    if (_dol_distrib == "exp") {
        float dol_mean = _disease.get_dol_mean();
        std::exponential_distribution<float> d(1.0/dol_mean);
        _dol_drawn = d(_RANDOM_GENERATOR);
        //		cout << "DEBUG: dol_drawn = " << _dol_drawn <<endl;
    }
    
    if (_doi_distrib == "dirac") _doi_drawn = _disease.get_doi_mean();
    if (_doi_distrib == "exp") {
        float doi_mean = _disease.get_doi_mean();
        std::exponential_distribution<float> d(1.0/doi_mean);
        _doi_drawn = d(_RANDOM_GENERATOR);
    }
    
    
    // TO DO: change that! Make this random...
    _is_symptomatic = true;
    _was_symptomatic = true;
}


void individual::recoverDisease() {
    /// This individual recovers from the disease
    
    _is_susceptible = false;
    _is_infected	= false;
    _is_infectious	= false;
    _is_recovered	= true;
    _is_symptomatic	= false;
    
    _doi= 0.0;
}


vector<individual> build_individuals(unsigned int n,
                                     const vector<schedule>& sched,
                                     string dol_distrib,
                                     string doi_distrib){
    /// Build several individuals
    /// TO DO: make it more sophisticated!
    
    vector<individual> x(n);
    
    std::uniform_real_distribution<double> unif(1.0, 80.0);
    std::uniform_real_distribution<double> unif_01(0.0, 1.0);
    std::uniform_int_distribution<unsigned long> unif_int(0.0, sched.size()-1);
    
    for (int i=0; i<n; i++)
    {
        double age = unif(_RANDOM_GENERATOR);
        individual tmp(i, age);
        
        tmp.set_immunity(0.0); //unif_01(_RANDOM_GENERATOR));
        tmp.set_frailty(1.0);  //unif_01(_RANDOM_GENERATOR));
        tmp.set_schedule(sched[unif_int(_RANDOM_GENERATOR)]);
        tmp.set_dol_distrib(dol_distrib);
        tmp.set_doi_distrib(doi_distrib);
        x[i] = tmp;
    }
    return x;
}


individual get_indiv_with_ID(ID id,
                             const vector<individual>& indiv_vec){
    /// Retrieve the first individual with ID='id' in a vector of individuals
    
    individual res;
    ID i=0;
    for (i=0; i<indiv_vec.size(); i++) {
        if (indiv_vec[i].get_id() == id) {
            res = indiv_vec[i];
            break;
        }
    }
    stopif(i>=indiv_vec.size(), "Individual ID " + to_string(id) + " not found in vector");
    
    return res;
}



ID individual::find_dest(unsigned int idx_timeslice){
    /// Find the ID of the social place this
    /// individual is supposed to go at a given times slice.
    
    SPtype sptype = get_schedule().get_sp_type()[idx_timeslice];
    
    //			Retrieve the actual destination:
    
    ID id_dest = __UNDEFINED_ID;
    if(sptype == SP_household)	id_dest = get_id_sp_household();
    if(sptype == SP_workplace)	id_dest = get_id_sp_workplace();
    if(sptype == SP_school)		id_dest = get_id_sp_school();
    if(sptype == SP_other)		id_dest = get_id_sp_other();
    if(sptype == SP_hospital)	id_dest = get_id_sp_hospital();
    if(sptype == SP_pubTransp)	id_dest = get_id_sp_pubTransp();
    
    return id_dest;
}









