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
	_isalive = true;
	_id = __UNDEFINED_ID;
	_age = 0.0;
	_immunity = 0.0;
	_frailty = 1.0;
	
	_id_sp_current		= __UNDEFINED_ID;
	_id_sp_household	= __UNDEFINED_ID;
	_id_sp_workplace	= __UNDEFINED_ID;
	_id_sp_school		= __UNDEFINED_ID;
	_id_sp_other		= __UNDEFINED_ID;
	_id_sp_hospital		= __UNDEFINED_ID;
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



void individual::set_id_sp_household(const socialPlace& sp){
	stopif(sp.get_type() != SP_household, "social space must be a household!");
	_id_sp_household = sp.get_id_sp();
}

void individual::set_id_sp_workplace(const socialPlace& sp){
	stopif(sp.get_type() != SP_workplace, "social space must be a workplace!");
	_id_sp_workplace = sp.get_id_sp();
}

void individual::set_id_sp_school(const socialPlace& sp){
	stopif(sp.get_type() != SP_school, "social space must be a school!");
	_id_sp_school = sp.get_id_sp();
}

void individual::set_id_sp_other(const socialPlace& sp){
	stopif(sp.get_type() != SP_other, "social space must be a other public space!");
	_id_sp_other = sp.get_id_sp();
}

void individual::set_id_sp_hospital(const socialPlace& sp){
	stopif(sp.get_type() != SP_hospital, "social space must be a hospital!");
	_id_sp_hospital = sp.get_id_sp();
}

void individual::set_id_sp_pubTransp(const socialPlace& sp){
	stopif(sp.get_type() != SP_pubTransp, "social space must be a public transportation!");
	_id_sp_pubTransp = sp.get_id_sp();
}



void individual::displayInfo(){
	
	cout << endl << "-- " << endl;
	cout << "individual ID: " << _id << endl;
	cout << "individual alive: " << _isalive << endl;
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
	cout << "-- " << endl;
}