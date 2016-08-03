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
    _id_sp_current      = __UNDEFINED_ID;
    _id_sp_household    = __UNDEFINED_ID;
    _id_sp_workplace    = __UNDEFINED_ID;
    _id_sp_school       = __UNDEFINED_ID;
    _id_sp_other        = __UNDEFINED_ID;
    _id_sp_hospital     = __UNDEFINED_ID;
    _id_sp_pubTransp    = __UNDEFINED_ID;
    
    // Clinical:
    _is_at_risk		= false;
    
    _is_susceptible     = true;
    _is_infected        = false;
    _is_latent          = false;
    _is_infectious      = false;
    _is_recovered       = false;
    
    _is_symptomatic     = false;
    _was_symptomatic    = false;

    _is_hosp            = false;
    _willbe_hosp        = false;
    _is_discharged      = false;
    _is_treated         = false;
    _is_vaccinated      = false;
    _will_die           = false;  
    
    
    _doi	= 0.0;
    _dol	= 0.0;
    _doh	= 0.0;
    _dobh   = 0.0;
    _doi_drawn	= 0.0;
    _dol_drawn	= 0.0;
    _doh_drawn	= 0.0;
    _dobh_drawn = 0.0;
    
    _dol_distrib = "dirac";
    _doi_distrib = "dirac";
    _doh_distrib = "dirac";
    
    _acquisition_time           = __UNDEFINED_FLOAT;
    _acquisition_time_infector  = __UNDEFINED_FLOAT;
    _num_secondary_cases        = 0;
    _ID_secondary_cases.clear();
}


individual::individual(){
    base_constructor();
}

individual::individual(ID id){
    base_constructor();
    _id = id;
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
//    	if(_id == 22){
//    		cout << "debug"<<endl;
//    	}
    // ------
    
    // Disease-related durations:
    if (_dol>0 && _dol <= _dol_drawn) _dol += dt;
    if (_doi>0 && _doi <= _doi_drawn) {
        _doi += dt;
        _dobh = _dol+_doi;
    }
    if (_doh>0 && _doh <= _doh_drawn) _doh += dt;
    
    // Update disease stages
    
    // From 'E' to 'I'
    if (_dol > _dol_drawn &&
        !_is_infectious &&
        !_is_recovered  &&
        _is_alive)
    {
        _is_latent      = false;
        _is_infectious  = true;
        _doi            = SUPERTINY;
        res = "E_to_I";
    }
    
    // form 'Is' to 'H'
    // * * * WARNING * * *
    // DO NOT DO ANYTHING FOR THIS CASE:
    // because hospitalization requires
    // moving the individual to an
    // 'hospital' social place,
    // so that's dealt with at the 'Simulation' level.
    
    
    // From 'I' to 'R'
    if (_doi > _doi_drawn &&
        !_is_recovered &&
        !_is_hosp      &&
        _is_alive)
    {
        // TO DO: use 'recoverDisease()' function instead
        
        _is_infectious  = false;
        _is_infected    = false;
        _is_symptomatic = false;
        _is_recovered   = true;
        res = "I_to_R";
    }
    
    // from 'H' to 'R'
    if (_doh > _doh_drawn &&
        !_is_discharged   &&
        _is_alive)
    {
        if (_will_die){
            //die();
            res = "DEATH";
        }
        else {
            _is_discharged = true;
            res = "DISCHARGED_HOSPITAL";
        }
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
    }
    
    if (_doi_distrib == "dirac") _doi_drawn = _disease.get_doi_mean();
    if (_doi_distrib == "exp") {
        float doi_mean = _disease.get_doi_mean();
        std::exponential_distribution<float> d(1.0/doi_mean);
        _doi_drawn = d(_RANDOM_GENERATOR);
    }
}

void individual::futureHospitalization(){
    /// This individual will be hopitalized,
    /// determine when and duration.
    /// Also determine if death occur at
    /// the end of the hospitalization period.
    
    _willbe_hosp = true;
    
    // hospitalization will happen after (drawn) latent,
    // but before end of (drawn) infectiousness
    std::uniform_real_distribution<> unif(0, _doi_drawn);
    _dobh_drawn = _dol_drawn + unif(_RANDOM_GENERATOR);
    
    float doh_mean = _disease.get_doh_mean();
    if (_doh_distrib == "dirac") _doh_drawn = doh_mean;
    if (_doh_distrib == "exp") {
        std::exponential_distribution<float> d(1.0/doh_mean);
        _doh_drawn = d(_RANDOM_GENERATOR);
    }
    // DEBUG
//    cout << "  DEBUG: doh_drawn = " << _doh_drawn << " " << _doh_distrib <<" mean: " <<doh_mean << endl;
}


void individual::recoverDisease() {
    /// This individual recovers from the disease
    
    _is_susceptible = false;
    _is_infected    = false;
    _is_infectious  = false;
    _is_recovered   = true;
    _is_symptomatic = false;
    
    _doi= 0.0;
}

void individual::futureDeath(){
    /// This individual will die at
    /// the end of the hospitalization period
    
    _will_die = true;
    
}

void individual::die(){
    /// This individual dies from the disease.
    /// NOTE: current implementation means death
    /// is only possible after hospitalization.
    
    stopif(!_is_hosp, "Individual cannot die if it's not already hospitalized!");
    
    _is_alive       = false;
    
    _is_hosp        = false;
    _is_susceptible = false;
    _is_infected    = false;
    _is_infectious  = false;
    _is_recovered   = false;
    _is_symptomatic = false;
    
    _doi = 0.0;
    _doh = 0.0;
}


vector<individual> build_individuals(uint n,
                                     const vector<schedule>& sched,
                                     string dol_distrib,
                                     string doi_distrib,
                                     string doh_distrib){
    /// Build several individuals
    /// TO DO: make it more sophisticated!
    
    vector<individual> x(n);
    
    std::uniform_real_distribution<double> unif_age(1.0, 80.0);
    std::uniform_real_distribution<double> unif_01(0.0, 1.0);
    std::uniform_int_distribution<unsigned long> unif_int(0.0, sched.size()-1);
    
    for (int i=0; i<n; i++)
    {
        double age = unif_age(_RANDOM_GENERATOR);
        individual tmp(i, age);
        
        tmp.set_immunity(unif_01(_RANDOM_GENERATOR)); //unif_01(_RANDOM_GENERATOR)
        tmp.set_frailty(unif_01(_RANDOM_GENERATOR));  //unif_01(_RANDOM_GENERATOR)
        tmp.set_schedule(sched[unif_int(_RANDOM_GENERATOR)]);
        tmp.set_dol_distrib(dol_distrib);
        tmp.set_doi_distrib(doi_distrib);
        tmp.set_doh_distrib(doh_distrib);
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



ID individual::find_dest(uint idx_timeslice){
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


void individual::reduce_doi_drawn(double x)
{
    _doi_drawn = _doi_drawn - x;
    if(_doi_drawn<0)
        _doi_drawn=SUPERTINY;
}


void individual::receive_treatment(double doi_reduction){
    /// This individual receives a treamtment
    /// that is supposed to reduced its DOI
    
    set_is_treated(true);
    std::exponential_distribution<double> expdist(1.0/doi_reduction);
    double tmp = _doi_drawn - expdist(_RANDOM_GENERATOR);
    double new_doi_drawn = (tmp>0)?tmp:SUPERTINY;
    set_doi_drawn(new_doi_drawn);

}


void individual::receive_vaccine(float time,
                                 float imm_incr,
                                 float frail_incr,
                                 float vaxlag){
    /// Vaccinate this individual
    
    _is_vaccinated = true;
    
    _imm_when_recv_vax      = _immunity;
    _frail_when_recv_vax    = _frailty;
    
    // TO DO: more sophisticated (age dependence?), investigate literature
    std::uniform_real_distribution<> unif_lag(vaxlag-vaxlag/2.0, vaxlag+vaxlag/2.0);
    std::uniform_real_distribution<> unif_imm(imm_incr-imm_incr/2.0, imm_incr+imm_incr/2.0);
    std::uniform_real_distribution<> unif_fra(frail_incr-frail_incr/2.0, frail_incr+frail_incr/2.0);
    
    _vax_time_received      = time;
    _vax_target_immunity    = _immunity + unif_imm(_RANDOM_GENERATOR);
    _vax_target_frailty     = _frailty + unif_fra(_RANDOM_GENERATOR);
    _vax_lag_full_efficacy  = unif_lag(_RANDOM_GENERATOR);
    
    if(_vax_target_immunity > 1.0) _vax_target_immunity = 1.0;
    if(_vax_target_frailty > 1.0)  _vax_target_frailty  = 1.0;
}


void individual::receive_cure(){
    /// Instantaneously cure this individual (--> doi_drawn=0)
    /// * * Used for debuging * *
    set_is_treated(true);
    set_doi_drawn(SUPERTINY);
}










