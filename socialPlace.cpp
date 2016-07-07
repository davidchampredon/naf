//
//  socialPlace.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-06.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "socialPlace.h"



string SPtype2string(SPtype x){
	
	string res = "";
	
	if(x == SP_household)	res = "Household";
	if(x == SP_school)		res = "School";
	if(x == SP_hospital)	res = "Hospital";
	if(x == SP_workplace)	res = "Workplace";
	if(x == SP_other)		res = "Other public space";
	if(x == SP_pubTransp)	res = "Public transport";
	if(x == SP_random)		res = "RANDOM";
	return res;
}


void socialPlace::base_constructor(){
	_size = 0;
	_prevalence = 0;
	_id_sp = __UNDEFINED_ID;
}

socialPlace::socialPlace(){
	base_constructor();
}

socialPlace::socialPlace(ID id, SPtype type){
	base_constructor();
	_id_sp = id;
	_type = type;
}

socialPlace::socialPlace(ID id_au, string name, ID id_region,
						 string regionName,
						 ID id_sp, SPtype type){
	base_constructor();
	// AU level:
	_id_au = id_au ;
	_name_au = name;
	_name_region = regionName;
	_id_region = id_region;
	
	// SP level:
	_id_sp = id_sp;
	_type = type;
}

socialPlace::socialPlace(areaUnit AU, ID id_sp, SPtype type){
	base_constructor();
	// AU level:
	_id_au			= AU.get_id_au();
	_name_au		= AU.get_name_au();
	_name_region	= AU.get_name_region();
	_id_region		= AU.get_id_region();
	
	// SP level:
	_id_sp = id_sp;
	_type = type;
}



void socialPlace::displayInfo(){
	cout << "--- " << endl ;
	cout << "SP type: " << _type << " (" << SPtype2string(_type) << ")" << endl;
	cout << "SP id: " << _id_sp << endl;
	cout << " Associated AU:";
	displayInfo_AU();

	cout << "Num. of indiv: " << _indiv.size() << endl;
	cout << "indiv ids: ";
	for(int i=0; i<_indiv.size(); i++) cout << _indiv[i].get_id() << "; ";
	cout << endl;
	cout << "SP prevalence: " << _prevalence <<endl;
	cout << "--- " << endl ;
}



void socialPlace::add_indiv(individual& newindiv){
	/// Add a new individual to _existing_ ones
	
	// update ID of this SP for new individual:
	newindiv.set_id_sp_current(_id_sp);
	
	_indiv.push_back(newindiv);
	// update size:
	_size++;
	
	// update prevalence:
	if(newindiv.is_infected()) _prevalence++;
}



void socialPlace::add_indiv(vector<individual>& newindiv){
	/// Add individuals to _existing_ ones

	for(int i=0; i<newindiv.size(); i++)
		add_indiv(newindiv[i]);
	
//	// set ID of this SP for all new individuals:
//	for(int i=0; i<newindiv.size(); i++) newindiv[i].set_id_sp_current(_id_sp);
//	
//	_indiv.insert(_indiv.end(), newindiv.begin(), newindiv.end());
	// update size:
//	_size = _indiv.size();
}


void socialPlace::remove_indiv(individual& x){
	/// Remove an individual from this social place
	
	// find the position of individual 'x' in the vector:
	auto pos = distance(_indiv.begin(), find(_indiv.begin(),_indiv.end(),x));
	stopif(pos>=_indiv.size(), "Try to remove ABSENT individual from social place!");

	_indiv.erase(_indiv.begin()+pos);

	// update ID (function 'add_individual' will specify ID)
	x.set_id_sp_current(__UNDEFINED_ID);
	
	// update size:
	_size--;
	// update prevalence:
	if(x.is_infected()) _prevalence--;
}


void socialPlace::remove_indiv(unsigned int pos){
	/// Remove an individual given its POSITION in vector '_indiv' in this social place.
	stopif(pos>=_indiv.size(), "Try to remove NON EXISTENT individual from social place!");
	
	bool is_infected = _indiv[pos].is_infected();
	
	_indiv.erase(_indiv.begin()+pos);
	
	// update size:
	_size--;
	// update prevalence:
	if(is_infected) _prevalence--;
}



void socialPlace::remove_indiv(vector<unsigned int> posvec){
	/// Remove SEVERAL individual given their INITIAL POSITION in vector '_indiv' in this social place.

//	stopif(posvec[max_element(posvec.begin(),posvec.end())]>=_indiv.size(), "Try to remove NON EXISTENT individual from social place!");

	// posvec must be sorted
	sort(posvec.begin(), posvec.end());
	
	for(int i=0; i<posvec.size(); i++){
		// update prevalence (must be before removing, else infection info is gone!):
		if(_indiv[posvec[i]-i].is_infected()) _prevalence--;
		// remove
		_indiv.erase(_indiv.begin()+posvec[i]-i); // '-i' to take into account the shrinking vector
		// update size:
		_size--;
	}
}




vector<unsigned int> socialPlace::pick_rnd_susceptibles(unsigned int num){
	
	/// Pick randomly 'num' susceptibles from this social place.
	/// Returns the POSITION of susceptibles in '_indiv' vector.
	
	stopif(_indiv.size() < num, "Asking for too many susceptibles!");
	
	// Shuffle elements (guarantees random pick)
	vector<individual> tmp = _indiv;
	random_shuffle(tmp.begin(), tmp.end());
	
	// Find the first 'num' susceptibles
	vector<unsigned int> pos;
	for (unsigned int i=0; pos.size() < num && i<tmp.size() ; i++) {
		if (tmp[i].is_infected()) pos.push_back(i);
	}
	
	return pos;	
}


unsigned int socialPlace::census_alive(){
	/// Counts all individuals alive
	
	unsigned int cnt = 0;
	for(int i=0; i<_indiv.size(); i++){
		if(_indiv[i].is_alive()) cnt++;
	}
	return cnt;
}
