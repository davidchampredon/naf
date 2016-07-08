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
	if(x == SP_RANDOM)		res = "RANDOM";
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
	
	cout << "Num. of indiv: " << _indiv.size() << " (check: " << _size << ")" << endl;
	cout << "indiv ids: ";
	for(int i=0; i<_indiv.size(); i++) cout << " [" <<_indiv[i].get_id() << "]";
	cout << endl;
	cout << "SP prevalence: " << _prevalence << " (check: "<< id_infected_bruteforce().size()  <<")";
	if(_prevalence != id_infected_bruteforce().size()) cout << " CHECK FAILED!";
	cout << endl;
	cout << "infected indiv ids: ";
	for(int i=0; i<_indiv.size(); i++)
		if(_indiv[i].is_infected()) {cout <<" ["<<_indiv[i].get_id() << "]";}
	
	cout << endl;
	cout << "--- " << endl ;
}



void socialPlace::add_indiv(individual & newindiv){
	/// Add a new individual to _existing_ ones
	
	// update ID of this SP for new individual:
	newindiv.set_id_sp_current(_id_sp);
	
	_indiv.push_back(newindiv);
	_size++;
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
	
	bool is_infected = _indiv[pos].is_infected(); // <-- MUST BE BEFORE erasing individual!
	
	_indiv.erase(_indiv.begin()+pos);
	
	// Mandatory updates:
	_size--;
	if(is_infected) _prevalence--;
}



void socialPlace::remove_indiv(vector<unsigned int> posvec){
	/// Remove SEVERAL individuals given their INITIAL POSITION in vector '_indiv' in this social place.
	
	//	stopif(posvec[max_element(posvec.begin(),posvec.end())]>=_indiv.size(), "Try to remove NON EXISTENT individual from social place!");
	
	// posvec must be sorted
	sort(posvec.begin(), posvec.end());
	
	for(int i=0; i<posvec.size(); i++){
		
		unsigned int idx = posvec[i]-i; // '-i' to take into account the shrinking vector
		
		cout << "DEBUG: removing indiv ID_"<<_indiv[idx].get_id();
		cout << " pos_"<<i<<" infected:"<<_indiv[idx].is_infected();
		
		// update prevalence (must be before removing, else infection info is gone!):
		if(_indiv[idx].is_infected()) _prevalence--;
		// remove
		_indiv.erase(_indiv.begin()+ idx );
		// update size:
		_size--;
		
		
		
		cout << " updated size: "<<_size<<" prev: "<<_prevalence << endl;
	}
}




vector<unsigned int> socialPlace::pick_rnd_susceptibles(unsigned int num){
	
	/// Pick randomly 'num' susceptibles from this social place.
	/// Returns the POSITION of susceptibles in '_indiv' vector.
	
	stopif(_indiv.size() < num, "Asking for too many susceptibles!");
	
	// WARNING: BRUTE FORCE => SLOW!
	// TO DO: OPTIMIZE
	
	vector<unsigned int> pos;
	
	for (unsigned int i=0; i<_indiv.size() ; i++){
		bool is_susceptible = !(_indiv[i].is_infected());
		if (is_susceptible) pos.push_back(i);
	}
	// Shuffle elements (guarantees random pick)
	random_shuffle(pos.begin(), pos.end());
	pos.resize(num);
	
	return pos;
	
	/* TRY TO OPTIMIZE BUT DOES NOT WORK!!!
	 
	 // Shuffle elements (guarantees random pick)
	 vector<individual> tmp = _indiv;
	 random_shuffle(tmp.begin(), tmp.end());
	 // Find the first 'num' susceptibles
	 vector<unsigned int> pos;
	 unsigned int cnt = 0;
	 for (unsigned int i=0; cnt < num && i<tmp.size() ; i++){
		bool is_susceptible = !(tmp[i].is_infected());
		if (is_susceptible) {
	 pos.push_back(i);
	 cnt++;
		}
	 }
	 return pos;
	 */
}


unsigned int socialPlace::census_alive(){
	/// Counts all individuals alive
	
	unsigned int cnt = 0;
	for(int i=0; i<_indiv.size(); i++){
		if(_indiv[i].is_alive()) cnt++;
	}
	return cnt;
}


vector<ID>	socialPlace::id_infected_bruteforce(){
	vector<ID> res;
	for(int i=0; i<_indiv.size(); i++){
		if(_indiv[i].is_infected()) res.push_back(_indiv[i].get_id());
	}
	return res;
}




ID socialPlace::find_dest(unsigned int pos, unsigned int idx_timeslice){
	/// Find the ID of the social place the individual is supposed to move to
	/// at the timeslice 'idx_timeslice' of the schedule.
	/// (individual is in position 'pos' in the vector '_indiv')
	// Retrieve the social place this individual is supposed
	// to go for this timeslice, according to its schedule
	SPtype sptype = _indiv[pos].get_schedule().get_sp_type()[idx_timeslice];
	
	//			Retrieve the actual destination:
	
	ID id_dest = __UNDEFINED_ID;
	if(sptype == SP_household)	id_dest = _indiv[pos].get_id_sp_household();
	if(sptype == SP_workplace)	id_dest = _indiv[pos].get_id_sp_workplace();
	if(sptype == SP_school)		id_dest = _indiv[pos].get_id_sp_school();
	if(sptype == SP_other)		id_dest = _indiv[pos].get_id_sp_other();
	if(sptype == SP_hospital)	id_dest = _indiv[pos].get_id_sp_hospital();
	if(sptype == SP_pubTransp)	id_dest = _indiv[pos].get_id_sp_pubTransp();
	
	return id_dest;
}


void socialPlace::acquireDisease(unsigned int pos){
	/// Individual at position 'pos' in '_indiv' acquires the disease
	
	_indiv[pos].acquireDisease();
	_prevalence++;
}



vector<socialPlace> build_random(unsigned int N, vector<areaUnit> auvec){
	/// Build randomly 'N' social places using provided area units
	
	unsigned long n_au = auvec.size();
	
	unsigned int nSPtype = SP_MAX;
	
	vector<socialPlace> v;
	
	for (int i=0; i<N; i++) {
		
		areaUnit A = auvec[rand()%n_au];				// choose randomly an AU
		SPtype sptype = (SPtype)(rand()%SP_MAX);	// choose randomly the type of SP
		socialPlace tmp(A,i,sptype);
		v.push_back(tmp);
		
		//DEBUG
		cout<<"DEBUG sptype = "<<SPtype2string(sptype)<< " ; area unit ID: " << A.get_id_au()<<endl;
	}
	
	return v;
}

void populate_random_with_indiv(vector<socialPlace> & sp,
								ID n_indiv,
								vector<schedule> sched){
	/// Create and distribute N individuals among all SP in 'sp'.
	/// Individuals have schedule randomly allocated.
	
	
	// STEP 1 - create individuals
	
	vector<individual> indivvec;
	
	for(int i=0; i<n_indiv; i++){
		double age = rand() % 90 + 1;
		
		individual tmp(i, age);
		
		unsigned int id_rnd;

		id_rnd = choose_SPtype_random(sp, SP_school);
		tmp.set_id_sp_school(sp[id_rnd]);
		
		id_rnd = choose_SPtype_random(sp, SP_hospital);
		tmp.set_id_sp_hospital(sp[id_rnd]);
		
		id_rnd = choose_SPtype_random(sp, SP_household);
		tmp.set_id_sp_household(sp[id_rnd]);
		
		id_rnd = choose_SPtype_random(sp, SP_workplace);
		tmp.set_id_sp_workplace(sp[id_rnd]);
		
		id_rnd = choose_SPtype_random(sp, SP_pubTransp);
		tmp.set_id_sp_pubTransp(sp[id_rnd]);
		
		id_rnd = choose_SPtype_random(sp, SP_other);
		tmp.set_id_sp_other(sp[id_rnd]);
		
		
		tmp.set_immunity(0.0);
		tmp.set_frailty(1.0);
		
		tmp.set_schedule(sched[rand() % sched.size()]);

		indivvec.push_back(tmp);
		// DEBUG
		// tmp.displayInfo();
	}
	
	// STEP 2 - assign individuals to a random SP
	
		for (int i=0; i<indivvec.size(); i++) {
			int sp_idx = rand() % sp.size();
			sp[sp_idx].add_indiv(indivvec[i]);
			//DEBUG
			// indivvec[i].displayInfo();
		}

}

unsigned int choose_SPtype_random(const vector<socialPlace>& sp, SPtype x){
	/// Choose randomly a SP with a given SPtype. Returns the position in the vector 'sp'
	
	vector<unsigned int > pos;
	for(unsigned int  i=0;i<sp.size(); i++){
		if(sp[i].get_type()==x) pos.push_back(i);
	}
	
	unsigned int choose_rnd = rand() % pos.size();
	unsigned int  res = pos[choose_rnd];
	return res;
}



void displayPopulationSize(const vector<socialPlace>& sp){
	
	cout<<endl<< " WORLD POPULATION "<<endl;
	ID s = 0;
	for (int i=0; i<sp.size(); i++) {
		ID sizei = sp[i].get_size();
		s+=sizei;
		cout << "pos= "<< i << " ; id= "<<sp[i].get_id_sp();
		cout << " ; popsize= "<< sizei;
		cout << " ; type= " << SPtype2string(sp[i].get_type()) <<endl;
	}
	cout << "Total population: "<< s << endl;
}









