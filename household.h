//
//  household.h
//  naf
//
//  Created by David CHAMPREDON on 2016-07-04.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#ifndef __naf__household__
#define __naf__household__

#include <stdio.h>
#include <vector>

#include "individual.h"


class household{
	
public:
	
	ID                 _id;
	unsigned long      _size;
	vector<ID>         _id_indiv;
	location           _location;
	
	household();
	household(ID id, location loc);
	
	
	location get_location(){return _location;}
	
	
	void populate_household(vector<individual>& indiv, vector<ID> idvec);
	
	void displayInfo();
	
};



#endif /* defined(__naf__household__) */
