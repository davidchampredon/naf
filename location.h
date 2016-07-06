//
//  location.h
//  naf
//
//  Created by David CHAMPREDON on 2016-07-04.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#ifndef __naf__location__
#define __naf__location__

#include <stdio.h>
#include <string>
#include "utils.h"

using namespace std;

class location{
	
public:
	string _name;
	ID     _id;
	
	location();
	location(ID id, string name) {
		_id = id ;
		_name = name;
	}
	
	void displayInfo();
	
};


#endif /* defined(__naf__location__) */
