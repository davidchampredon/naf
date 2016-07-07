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
	bool	_is_infected;
	double	_doi; // duration of infection
	
	// Social spaces linked to this individual
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
	
	void set_id_sp_household(const socialPlace& sp);
	void set_id_sp_workplace(const socialPlace& sp);
	void set_id_sp_school(const socialPlace& sp);
	void set_id_sp_other(const socialPlace& sp);
	void set_id_sp_hospital(const socialPlace& sp);
	void set_id_sp_pubTransp(const socialPlace& sp);

	void set_immunity(double x)	{_immunity = x;}
	void set_frailty(double x)	{_frailty = x;}
	
	void set_schedule(schedule s) {_schedule = s;}
	
	void forget_id_sp_household(){_id_sp_household = __UNDEFINED_ID;}
	
	
	// Get functions
	
	ID get_id(){return _id;}
	ID get_id_sp_household()	{return _id_sp_household;}
	ID get_id_sp_workplace()	{return _id_sp_workplace;}
	ID get_id_sp_school()		{return _id_sp_school;}
	ID get_id_sp_other()		{return _id_sp_other;}
	ID get_id_sp_hospital()		{return _id_sp_hospital;}
	ID get_id_sp_pubTransp()	{return _id_sp_pubTransp;}
	
	double	get_immunity()		{return _immunity;}
	double	get_frailty()		{return _frailty;}
	
	bool is_infected()		{return _is_infected;}
	bool is_alive()			{return _is_alive;}
	
	schedule get_schedule() {return _schedule;}
	
	
	// Epidemiology
	void acquireDisease() {_is_infected = true; _doi= SUPERTINY;}
	void recoverDisease() {_is_infected = false; _doi= 0.0;}
	
	// Miscellenaous
	void displayInfo();
};

// Operators
inline bool operator == ( individual a, individual b){
	return (a.get_id() == b.get_id());
}



#endif /* defined(__naf__individual__) */
