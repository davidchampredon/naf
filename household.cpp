//
//  household.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-04.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "household.h"



household::household(){
	_id = __UNDEFINED_ID;
	_size = 0;
	
	location none;
	_location = none;
}


household::household(ID id, location loc){
	_id = id;
	_location = loc;
}



void household::populate_household(vector<individual>& indiv, vector<ID> idvec){
	
	// link hh to indiv
	_id_indiv.clear();
	for(ID i=0; i<idvec.size(); i++) _id_indiv.push_back(indiv[idvec[i]]._id);
	_size = idvec.size();
	
	// link indiv to hh
	for(ID i=0; i<idvec.size(); i++) indiv[idvec[i]]._id_household = _id;
}


void household::displayInfo(){
	cout << endl;
	cout << "Household ID: " << _id << endl;
	cout << "Household size: " << _size << endl;
	cout << "Individuals ID: ";
	for(int i=0; i<_id_indiv.size(); i++) cout<<_id_indiv[i]<<" ; ";
//	cout << endl;
	_location.displayInfo();
}