//
//  areaUnit.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-06.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "areaUnit.h"


#include <iostream>

areaUnit::areaUnit(){
	_id_au = __UNDEFINED_ID;
	_name_au = "UNASSIGNED_AREAUNIT_NAME";
}



void areaUnit::displayInfo(){
	cout << endl;
	cout << "areaUnit ID:   " << _id_au   << endl;
	cout << "areaUnit name: " << _name_au << endl;
	cout << "areaUnit's region: ID="<<_id_region<<" ; name="<<_name_region<<endl;
}

void areaUnit::displayInfo_AU(){
	displayInfo();
}