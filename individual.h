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
#include "location.h"
#include "utils.h"


class individual{

public:
	ID     _id;
	
	// Biological
	double _age;
	double _immunity;   // between 0 and 1.0
	double _frailty;	// between 0 and 1.0
	
	// Spatial
	location _current_location;
	ID       _id_household;
	
	location _work_location;
	location _school_location;
	
	// Constructors
	individual();
	individual(ID id, double age);
	individual(ID id, double age, ID id_household);
	
	// Set functions
	void set_id_household(ID id_hh){_id_household = id_hh;}
	
	void displayInfo();
};



#endif /* defined(__naf__individual__) */
