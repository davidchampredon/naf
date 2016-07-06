//
//  individual.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-04.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "individual.h"



individual::individual(){
	
	_id = __UNDEFINED_ID;
	_age = 0.0;
	_immunity = 0.0;
	_frailty = 1.0;
}


individual::individual(ID id, double age){
	_id = id;
	_age = age;
}

individual::individual(ID id, double age, ID id_household){
	_id = id;
	_age = age;
	_id_household = id_household;
}


void individual::displayInfo(){
	
	cout << endl << " ---- Individual Info ---- " << endl;
	cout << "individual ID: " << _id << endl;
	cout << "individual age: " << _age << endl;
	cout << "individual's current location: " << _current_location._name << endl;
	cout << "individual's household id: " << _id_household << endl;
	
}