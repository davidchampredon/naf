//
//  location.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-04.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "location.h"
#include <iostream>

location::location(){
	_id = __UNDEFINED_ID;
	_name = "UNASSIGNED_LOCATION_NAME";
}



void location::displayInfo(){
	cout << endl;
	cout << "Location ID:   " << _id   << endl;
	cout << "Location name: " << _name << endl;
}