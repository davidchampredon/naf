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



string SPtype2string(SPtype x);


class socialPlace: public areaUnit{

protected:
	
	SPtype			_type;
	ID				_id_sp;
	unsigned long	_size;
	
	vector<individual> _indiv;

	unsigned int	_prevalence;
	
	
public:
	
	// Constructors:
	void base_constructor();
	socialPlace();
	socialPlace(ID id, SPtype type);
	socialPlace(ID id_au, string name, ID id_region, string regionName,
				ID id_sp, SPtype type);
	socialPlace(areaUnit AU, ID id_sp, SPtype type);
	
	
	// Operations on individuals:
	void add_indiv(individual& indiv);
	void add_indiv(vector<individual>& indiv);
	void remove_indiv(individual& indiv);
	void remove_indiv(unsigned int pos);
	void remove_indiv(vector<unsigned int> posvec);
	
	// Set functions
	void set_prevalence(unsigned int p) {_prevalence = p;}
	void increase_prevalence() {_prevalence++;}
	
	// Get functions:
	SPtype	get_type() const {return _type;}
	ID		get_id_sp() const {return _id_sp;}
	unsigned long			get_size(){return _size;}
	unsigned int			get_prevalence(){return _prevalence;}
	vector<individual>	get_indiv() {return _indiv;}
	
	
	// Diseases:
	vector<unsigned int> pick_rnd_susceptibles(unsigned int num);
	
	
	// Miscellenaous:
	void displayInfo();
	unsigned int	census_alive();
	
};


#endif /* defined(__naf__socialPlace__) */
