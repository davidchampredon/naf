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

#include "utils.h"
#include "schedule.h"


enum scheduleType {
	workerSedentary, workerTravel,
	traveler,
	stayHome, hospitalized,
	leisure};

class socialPlace;

class individual{

protected:

	ID		_id;
	double	_age;
	bool	_is_alive;
	
	// Medical
	double	_immunity;   // between 0 and 1.0 (1=completely immune)
	double	_frailty;	// between 0 and 1.0 (1=extremely weak)
	
	// Disease
	bool	_is_latent;
	
	bool	_is_infected;
	bool	_is_infectious;
	bool	_is_symptomatic;

	bool	_is_recovered;
	bool	_is_hosp;
	
	double	_dol; // duration of latency
	double	_doi; // duration of infection
	double	_doh; // duration of hospitalization

	// Risk group:
	// has an underlying condition increasing
	// the risk of severity if infected by disease
	// (e.g. HIV, pregnant, etc.)
	
	// ???? REDUNDANT WITH FRAILTY ????
	bool	_is_at_risk;
	
	
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
	
	
public:
	
	// Constructors
	
	void base_constructor();
	individual();
	individual(ID id, double age);
	individual(ID id, double age, ID id_household);
	
	
	// Set functions
	
	void set_id_sp_current(ID id_sp){_id_sp_current = id_sp;}
	
	void set_id_sp_household(socialPlace& sp);
	void set_id_sp_workplace(socialPlace& sp);
	void set_id_sp_school(socialPlace& sp);
	void set_id_sp_other(socialPlace& sp);
	void set_id_sp_hospital(socialPlace& sp);
	void set_id_sp_pubTransp(socialPlace& sp);
	
	void set_id_sp(SPtype type, socialPlace& sp);

	void set_immunity(double x)	{_immunity = x;}
	void set_frailty(double x)	{_frailty = x;}
	
	void set_schedule(schedule s) {_schedule = s;}
	
	void forget_id_sp_household(){_id_sp_household = __UNDEFINED_ID;}
	void set_is_infected(bool x) {_is_infected = x;}
	
	
	// Get functions
	
	ID get_id()					const {return _id;}
	ID get_id_sp_household()	const {return _id_sp_household;}
	ID get_id_sp_workplace()	const {return _id_sp_workplace;}
	ID get_id_sp_school()		const {return _id_sp_school;}
	ID get_id_sp_other()		const {return _id_sp_other;}
	ID get_id_sp_hospital()		const {return _id_sp_hospital;}
	ID get_id_sp_pubTransp()	const {return _id_sp_pubTransp;}
	
	double	get_immunity()		const {return _immunity;}
	double	get_frailty()		const {return _frailty;}
	
	bool is_infected()			const {return _is_infected;}
	bool is_infectious()		const {return _is_infectious;}
	bool is_symptomatic()		const {return _is_symptomatic;}
	bool is_latent()			const {return _is_latent;}
	bool is_recovered()			const {return _is_recovered;}
	bool is_hosp()				const {return _is_hosp;}
	
	bool is_alive()				const {return _is_alive;}
	
	schedule get_schedule() {return _schedule;}
	
	
	// Epidemiology
	double calc_proba_acquire_disease();
	void acquireDisease();
	void recoverDisease() {_is_infected = false; _doi= 0.0;}
	
	// Miscellenaous
	ID find_dest(unsigned int idx_timeslice);
	void displayInfo();
};

// Operators
inline bool operator == ( individual a, individual b){
	return (a.get_id() == b.get_id());
}


inline void acquireDisease(individual& x) {x.set_is_infected(true);}

vector<individual> build_individuals(unsigned int n, const vector<schedule>& sched);


individual get_indiv_with_ID(ID id, const vector<individual>& indiv_vec);


#endif /* defined(__naf__individual__) */
