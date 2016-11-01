//
//  individual.h
//  naf
//
//  Created by David CHAMPREDON on 2016-07-04.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#ifndef __naf__individual__
#define __naf__individual__

#include <stdio.h>
#include <iostream>

#include "dcTools.h"
#include "schedule.h"
#include "disease.h"


enum scheduleType {
	workerSedentary,
	workerTravel,
	traveler,
	stayHome,
	hospitalized,
	leisure
};

class socialPlace;

class individual{

protected:

	ID		_id;
	float	_age;
	bool	_is_alive;
	
	// Medical
	float	_immunity_hum;  // humoral immunity index - between 0 and 1.0 (1=completely immune)
    float	_immunity_cell; // cellular immunity index - between 0 and 1.0 (1=completely immune)
	float	_frailty;		// between 0 and 1.0 (1=extremely weak)
	
	// Disease
	
	disease	_disease;
	
	bool    _is_susceptible;
	bool    _is_latent;
	bool    _is_infected;
	bool    _is_infectious;
	bool    _is_symptomatic;

    bool    _is_recovered;
    bool    _is_hosp;
    bool    _willbe_hosp; // hospitalization is decided in advance, at transmission time
    bool    _is_discharged;
    bool    _will_die;      // whether death occurs at the end of hospitalization.
    
    bool    _was_symptomatic;   // records if was symptomatic (when checked after infection).
    bool    _was_hosp;          // records if has been hospitalized
    
    float   _dol;  // duration of latency
    float   _doi;  // duration of infection
    float   _dobh; // duration before hospitalization
    float   _doh;  // duration of hospitalization
	
    // When individual is infected,
	// disease stages durations are randomly drawn:
    float	_dol_drawn;
    float	_doi_drawn;
    float	_doh_drawn;
    float	_dobh_drawn;
    
    string	_dol_distrib;  // distribution of the duration of latency
    string	_doi_distrib;  // distribution of the duration of infectiousness
    string	_doh_distrib;  // distribution of the duration of hospitalization
	
	// ???? REDUNDANT WITH FRAILTY ????
	bool	_is_at_risk;
	
	// Epidemiology at individual level
    float       _acquisition_time;
    float       _acquisition_time_infector;
    
    uint        _num_secondary_cases;
    vector<ID>  _ID_secondary_cases;

    // Treatmentand Vaccination
    bool    _is_treated;
    bool    _is_vaccinated;
    
    // After being vaccinated,
    // there is a lag until full
    // vaccine efficacy.
    // These variables represent the
    // time lag and the target levels
    // reached at the end of that period:
    float   _vax_time_received;
    float   _vax_lag_full_efficacy; // <-- lag b/w time vax recvd and full efficacy
    float   _vax_target_immunity_hum;
    float   _vax_target_immunity_cell;
    float   _vax_target_frailty;
    float   _imm_hum_when_recv_vax;      // humoral immunity when vaccinated
    float   _imm_cell_when_recv_vax;     // cellular immunity when vaccinated
    float   _frail_when_recv_vax;        // frailty when vaccinated
    
	// Social spaces linked
	// to this individual
	ID	_id_sp_current; // useful?
	ID	_id_sp_household;
	ID	_id_sp_workplace;
	ID	_id_sp_school;
	ID	_id_sp_other;
	ID	_id_sp_hospital;
	ID	_id_sp_pubTransp;
	
	// Schedule
	schedule _schedule;
	
	// other private functions:
	
	string	_disease_status_update(double dt);
	
	
public:
	
	// Constructors
	
	void base_constructor();
	individual();
     individual(ID id);
	individual(ID id, float age);
	individual(ID id, float age, ID id_household);
	
	
	// Set functions
	
    void set_id(uint x) {_id = x;}
    void set_age(double x){_age = x;}
    
	void set_id_sp_current(ID id_sp){_id_sp_current = id_sp;}
	
	void set_id_sp_household(socialPlace& sp);
	void set_id_sp_workplace(socialPlace& sp);
	void set_id_sp_school(socialPlace& sp);
	void set_id_sp_other(socialPlace& sp);
	void set_id_sp_hospital(socialPlace& sp);
	void set_id_sp_pubTransp(socialPlace& sp);
	void set_id_sp(SPtype type, socialPlace& sp);

	void forget_id_sp_household(){_id_sp_household = __UNDEFINED_ID;}
	
	void set_immunity_hum(float x)      {_immunity_hum = x;}
    void set_immunity_cell(float x)     {_immunity_cell = x;}
	void set_frailty(float x)           {_frailty = x;}
	
	void set_schedule(schedule s)       {_schedule = s;}
	
	void set_disease(const disease& d)  {_disease = d;}
	void set_is_infected(bool x)        {_is_infected = x;}
    void set_is_hosp(bool x)            {_is_hosp = x;}
    void set_willbe_hosp(bool x)        {_willbe_hosp = x;}
    void set_is_discharged(bool x)      {_is_discharged = x;}
    
    void set_is_symptomatic(bool x)     {_is_symptomatic = x;}
    void set_was_symptomatic(bool x)    {_was_symptomatic = x;}
    void set_was_hosp(bool x)           {_was_hosp = x;}
    
    void set_is_treated(bool x)         {_is_treated = x;}
    void set_is_vaccinated(bool x)      {_is_vaccinated = x;}
    
    void set_doi(double x)              {_doi = x;}
    void set_doi_drawn(double x)        {_doi_drawn = x;}
    void reduce_doi_drawn(double x);
    void set_doh(double x)              {_doh = x;}
    
    void set_dol_distrib(string d)      {_dol_distrib = d;}
    void set_doi_distrib(string d)      {_doi_distrib = d;}
    void set_doh_distrib(string d)      {_doh_distrib = d;}
    
    void set_acquisition_time(float t)  {_acquisition_time = t;}
    void set_acquisition_time_infector(float t)  {_acquisition_time_infector = t;}
    
    void set_num_secondary_cases(uint x)    {_num_secondary_cases = x;}
    void increment_num_secondary_cases()    {_num_secondary_cases++;}
    
    void set_ID_secondary_cases(vector<ID> x)   {_ID_secondary_cases = x;}
    void push_ID_secondary_cases(ID x)          {_ID_secondary_cases.push_back(x);}
    
    
	// Get functions
	
	ID get_id()                 const {return _id;}
	ID get_id_sp_household()    const {return _id_sp_household;}
	ID get_id_sp_workplace()    const {return _id_sp_workplace;}
	ID get_id_sp_school()       const {return _id_sp_school;}
	ID get_id_sp_other()        const {return _id_sp_other;}
	ID get_id_sp_hospital()     const {return _id_sp_hospital;}
	ID get_id_sp_pubTransp()    const {return _id_sp_pubTransp;}
	
    double	get_age()           const {return _age;}
    float	get_immunity_hum()  const {return _immunity_hum;}
    float	get_immunity_cell() const {return _immunity_cell;}
    float	get_frailty()       const {return _frailty;}
    float	get_dol_drawn()     const {return _dol_drawn;}
    float	get_doi_drawn()     const {return _doi_drawn;}
    float	get_doh_drawn()     const {return _doh_drawn;}
    float	get_dobh_drawn()    const {return _dobh_drawn;}
    float	get_dol()           const {return _dol;}
    float	get_doi()           const {return _doi;}
    float	get_doh()           const {return _doh;}
    float	get_dobh()          const {return _dobh;}
    
    float get_vax_time_received()       const {return _vax_time_received;}
    float get_vax_lag_full_efficacy()   const {return _vax_lag_full_efficacy;}
    float get_vax_target_immunity_hum() const {return _vax_target_immunity_hum;}
    float get_vax_target_immunity_cell()const {return _vax_target_immunity_cell;}
    float get_vax_target_frailty()      const {return _vax_target_frailty;}
    float get_imm_hum_when_recv_vax()   const {return _imm_hum_when_recv_vax;}
    float get_imm_cell_when_recv_vax()  const {return _imm_cell_when_recv_vax;}
    float get_frail_when_recv_vax()     const {return _frail_when_recv_vax;}
	
    float get_acquisition_time()            const{return _acquisition_time;}
    float get_acquisition_time_infector()   const{return _acquisition_time_infector;}
    uint  get_num_secondary_cases()         const{return _num_secondary_cases;}
    
    bool is_susceptible()   const {return _is_susceptible;}
    bool is_infected()      const {return _is_infected;}
    bool is_infectious()    const {return _is_infectious;}
    bool is_symptomatic()   const {return _is_symptomatic;}
    bool was_symptomatic()  const {return _was_symptomatic;}
    bool was_hosp()         const {return _was_hosp;}
    bool is_latent()        const {return _is_latent;}
    bool is_recovered()     const {return _is_recovered;}
    bool is_hosp()          const {return _is_hosp;}
    bool willbe_hosp()      const {return _willbe_hosp;}
    bool will_die()         const {return _will_die;}
    bool is_discharged()    const {return _is_discharged;}
    bool is_treated()       const {return _is_treated;}
    bool is_vaccinated()    const {return _is_vaccinated;}
	
	bool is_alive()				const {return _is_alive;}
	
	schedule get_schedule()     const {return _schedule;}
	
    disease get_disease()       const {return _disease;}
    
    SPtype get_schedule_sp_type(uint idx_timeslice) const{
        return _schedule.get_sp_type(idx_timeslice);
    }
	
	// Time updates for all relevant members.
	string	time_update(double dt);
    bool    time_to_hospitalize() {return (_dobh >= _dobh_drawn);}
	
	// Epidemiology
    
    void    die();
    void    acquireDisease();
    void    recoverDisease();
    void    futureHospitalization();
    void    futureDeath();
    void    receive_treatment(double doi_reduction);
    void    receive_vaccine(float time,
                            float vaxlag,
                            float imm_hum_incr,
                            float imm_cell_incr,
                            float frail_incr);
    void    receive_cure();
	
	// Miscellenaous
	ID find_dest(uint idx_timeslice, bool print_error = false);
	void displayInfo();
};

// Operators
inline bool operator == ( individual a, individual b){
	return (a.get_id() == b.get_id());
}


//DELETE WHEN SURE: inline void acquireDisease(individual& x) {x.set_is_infected(true);}

vector<individual> build_individuals(uint n,
									 const vector<schedule>& sched,
									 string dol_distrib,
									 string doi_distrib,
                                     string doh_distrib);


individual get_indiv_with_ID(ID id, const vector<individual>& indiv_vec);


#endif /* defined(__naf__individual__) */
