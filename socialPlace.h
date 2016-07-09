//
//  socialPlace.h
//  naf
//
//  Created by David CHAMPREDON on 2016-07-06.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#ifndef __naf__socialPlace__
#define __naf__socialPlace__

#include <stdio.h>


#include "areaUnit.h"
#include "individual.h"
#include "utils.h"
#include "probaDistribution.h"


string SPtype2string(SPtype x);
SPtype int2SPtype(unsigned int i);


class socialPlace: public areaUnit{

protected:
	
	SPtype			_type;
	ID				_id_sp;  
	unsigned long	_size;
	
	vector<individual> _indiv;
	
	vector<ID>		_linked_indiv_id;  // ID of individuals that are linked to this social place

	unsigned int	_prevalence;
	
	
public:
	
	// Constructors:
	void base_constructor();
	socialPlace();
	socialPlace(ID id, SPtype type);
	socialPlace(ID id_au, string name, ID id_region, string regionName,
				ID id_sp, SPtype type);
	socialPlace(areaUnit AU, ID id_sp, SPtype type);
	
	
	// Movements of individuals:
	ID find_dest(unsigned int pos, unsigned int idx_timeslice);
	void add_indiv(individual& indiv);
	void add_indiv(vector<individual>& indiv);
	void remove_indiv(individual& indiv);
	void remove_indiv(unsigned int pos);
	void remove_indiv(vector<unsigned int> posvec);
	void add_linked_indiv(ID id);
	void remove_linked_indiv(ID id);
	ID   n_linked_indiv(){return _linked_indiv_id.size();}
	
	// Set functions
	void set_prevalence(unsigned int p) {_prevalence = p;}
	void increase_prevalence() {_prevalence++;}
	void set_type(SPtype t){_type = t;}
	
	
	// Get functions:
	SPtype			get_type() const {return _type;}
	ID				get_id_sp() const {return _id_sp;}
	unsigned long	get_size() const {return _size;}
	unsigned int	get_prevalence() const {return _prevalence;}
	vector<individual>	get_indiv() {return _indiv;}
	individual		get_indiv(unsigned int pos) {return _indiv[pos];}
	vector<ID>		get_linked_indiv_id() {return _linked_indiv_id;}
	
	
	// Diseases:
	void acquireDisease(unsigned int pos);
	vector<unsigned int> pick_rnd_susceptibles(unsigned int num);
	vector<ID>	id_infected_bruteforce();
	
	
	// Miscellenaous:
	void displayInfo();
	unsigned int	census_alive();
};


vector<socialPlace> build_world_random(unsigned int n_sp, vector<areaUnit> auvec);
vector<socialPlace> build_world_simple(vector<SPtype> spt,
									   vector<unsigned int> n_sp,
									   vector< probaDistrib<unsigned int> > p_size,
									   vector<individual>& indiv,
									   vector<areaUnit> auvec,
									   unsigned int seed =12345);

void populate_random_with_indiv(vector<socialPlace>& v, unsigned int total_indiv, vector<schedule> sched);

unsigned int  choose_SPtype_random(const vector<socialPlace>& sp, SPtype x);


void displayPopulationSize(const vector<socialPlace>& sp);





#endif /* defined(__naf__socialPlace__) */
