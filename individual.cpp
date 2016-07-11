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
	_is_alive = true;
	_id = __UNDEFINED_ID;
	_age = 0.0;
	_immunity = 0.0;
	_frailty = 1.0;
	
	_id_sp_current		= __UNDEFINED_ID;
	_id_sp_household	= __UNDEFINED_ID;
	_id_sp_workplace	= __UNDEFINED_ID;
	_id_sp_school		= __UNDEFINED_ID;
	_id_sp_other		= __UNDEFINED_ID;
	_id_sp_hospital	= __UNDEFINED_ID;
	_id_sp_pubTransp	= __UNDEFINED_ID;
	
	_is_infected = false;
	_doi = 0.0;
}


individual::individual(){
	base_constructor();
}


individual::individual(ID id, double age){
	base_constructor();
	_id = id;
	_age = age;
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
	cout << "-- " << endl;
}




double individual::calc_proba_acquire_disease(){
	/// Calculate the probability of aquiring
	/// the disease, given an infectious contact
	
	// TO DO: maybe implement something more sophisticated ? (logistic curves?)
	
	return (1-_immunity) * _frailty ;
}

void individual::acquireDisease(){
	 _is_infected = true;
	_doi = SUPERTINY;
}



vector<individual> build_individuals(unsigned int n, const vector<schedule>& sched){
	/// Build several individuals
	/// TO DO: make it more sophisticated!
	
	vector<individual> x(n);
	
	for (int i=0; i<n; i++) {
		double age = rand() % 80 + 1;
		individual tmp(i, age);
		
		tmp.set_immunity((double) rand() / (RAND_MAX));
		tmp.set_frailty((double) rand() / (RAND_MAX));
		
		tmp.set_schedule(sched[rand() % sched.size()]);
		
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









