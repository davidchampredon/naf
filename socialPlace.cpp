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
	if(x == SP_school)	res = "School";
	if(x == SP_hospital)	res = "Hospital";
	if(x == SP_workplace)	res = "Workplace";
	if(x == SP_other)		res = "Other public space";
	if(x == SP_pubTransp)	res = "Public transport";
	if(x == SP_RANDOM)	res = "RANDOM";
	return res;
}


SPtype int2SPtype(unsigned int i){
	// warning: order matters!
	SPtype res = SP_MAX;
	
	if (i==0) res = SP_household;
	if (i==1) res = SP_school;
	if (i==2) res = SP_hospital;
	if (i==3) res = SP_workplace;
	if (i==4) res = SP_other;
	if (i==5) res = SP_pubTransp;
	
	return res;
}

void socialPlace::base_constructor(){
	_size = 0;
	
	_prevalence = 0;
	_n_S = _size;
	_n_E = 0;
	_n_Ia = 0;
	_n_Is = 0;
	_n_D = 0;
	_n_H = 0;
	_n_R = 0;
	
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
	_id_sp = id_sp;
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



void socialPlace::update_epidemic_count(const individual& indiv,
										string update_type){
	/// Update epidemic count when a given
	/// individual is either added or removed
	
	
	// Updates following a movement
	// from one social place to another
	
	short change = 0;
	
	if (update_type == "add")		change = +1;
	if (update_type == "remove")	change = -1;
	
	if (change != 0) {
		if ( indiv.is_infected() )
		{
			_prevalence = _prevalence + change;
			
			if( indiv.is_latent()) {
				_n_E = _n_E + change;
			}
			
			if( indiv.is_infectious() )
			{
				if( indiv.is_symptomatic() )
				{ _n_Is = _n_Is + change; }
				else
				{ _n_Ia = _n_Ia + change; }
			}
		}
		
		if (!indiv.is_infected()){
			if (indiv.is_recovered())	_n_R += change;
			if (!indiv.is_recovered())	_n_S += change;
		}
	}
	
	// Updates following new case infected.
	// This MUST be AFTER the 'if' above:
	if (update_type == "new_case"){
		_n_E += 1;
		_n_S += -1;
	}
	
}


void socialPlace::add_indiv(individual & newindiv){
	/// Add a new individual to _existing_ ones
	
	// update ID of this SP for new individual:
	newindiv.set_id_sp_current(_id_sp);
	
	_indiv.push_back(newindiv);
	_size++;
	// DELETE WHEN SURE: if(newindiv.is_infected()) _prevalence++;
	update_epidemic_count(newindiv, "add");
}


void socialPlace::add_indiv(vector<individual>& newindiv){
	/// Add individuals to _existing_ ones
	
	for(int i=0; i<newindiv.size(); i++)
		add_indiv(newindiv[i]);
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
	// DELETE WHEN SURE: if(x.is_infected()) _prevalence--;
	update_epidemic_count(x, "remove");
}


void socialPlace::remove_indiv(unsigned int pos){
	/// Remove an individual given its POSITION
	/// in vector '_indiv' in this social place.
	
	stopif(pos>=_indiv.size(), "Try to remove NON EXISTENT individual from social place!");
	
	// DELETE WHEN SURE: bool is_infected = _indiv[pos].is_infected(); // <-- MUST BE BEFORE erasing individual!
	update_epidemic_count(_indiv[pos], "remove"); // <-- MUST BE BEFORE erasing individual!
	
	_indiv.erase(_indiv.begin()+pos);
	
	// Mandatory updates:
	_size--;
	// DELETE WHEN SURE: if(is_infected) _prevalence--;
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
		// DELETE WHEN SURE: if(_indiv[idx].is_infected()) _prevalence--;
		update_epidemic_count(_indiv[idx], "remove");
		
		// remove
		_indiv.erase(_indiv.begin()+ idx );
		// update size:
		_size--;

		cout << " updated size: "<<_size<<" prev: "<<_prevalence << endl;
	}
}

void socialPlace::add_linked_indiv(ID id){
	/// Add the ID of an individual who is supposed to be linked to this social place
	_linked_indiv_id.push_back(id);
}

void socialPlace::remove_linked_indiv(ID id){
	/// Remove the ID of an individual who is supposed to be linked to this social place
	_linked_indiv_id.erase(std::remove(_linked_indiv_id.begin(), _linked_indiv_id.end(), id));
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


unsigned int socialPlace::census_infectious(){
	/// Counts all infectious individuals (brute force, hence slow!)
	
	unsigned int cnt = 0;
	for(int i=0; i<_indiv.size(); i++){
		if(_indiv[i].is_infectious()) cnt++;
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
	
	return _indiv[pos].find_dest(idx_timeslice);
}



ID socialPlace::find_dest_linked(unsigned int pos,
								 unsigned int idx_timeslice,
								 const vector<individual>& indiv_vec){
	/// Find the ID of the social place the LINKED individual is supposed to move to
	/// at the timeslice 'idx_timeslice' of the schedule.
	/// (individual is in position 'pos' in the vector of LINKED individuals)
	
	individual tmp = get_indiv_with_ID(_linked_indiv_id[pos], indiv_vec);
	
	return tmp.find_dest(idx_timeslice);
}



void socialPlace::set_disease_to_all_indiv(const disease & d){
	/// Set the disease 'd' to all individuals in the social place
	
	for (ID i=0; i<_indiv.size(); i++) {
		_indiv[i].set_disease(d);
	}
}


void socialPlace::acquireDisease(unsigned int pos){
	/// Individual at position 'pos' in '_indiv' acquires the disease
	
	_indiv[pos].acquireDisease();
	update_epidemic_count(_indiv[pos], "new_case");
}



vector<socialPlace> build_world_random(unsigned int N, vector<areaUnit> auvec){
	/// Build randomly 'N' social places using provided area units
	
	unsigned long n_au = auvec.size();
	vector<socialPlace> v;
	
	for (int i=0; i<N; i++) {
		
		areaUnit A = auvec[rand()%n_au];				// choose randomly an AU
		SPtype sptype = (SPtype)(rand()%SP_MAX);	// choose randomly the type of SP
		socialPlace tmp(A,i,sptype);
		v.push_back(tmp);
		
		//DEBUG
		//cout<<"DEBUG sptype = "<<SPtype2string(sptype)<< " ; area unit ID: " << A.get_id_au()<<endl;
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
	unsigned int res = pos[choose_rnd];
	return res;
}



void displayPopulationSize(const vector<socialPlace>& sp){
	
	cout<<endl<< " WORLD POPULATION "<<endl;
	ID s = 0;
	for (int i=0; i<sp.size(); i++) {
		ID sizei = (ID)(sp[i].get_size());
		s+=sizei;
		cout << "pos= "<< i << " ; id= "<<sp[i].get_id_sp();
		cout << " ; popsize= "<< sizei;
		cout << " ; type= " << SPtype2string(sp[i].get_type()) <<endl;
	}
	cout << "Total population: "<< s << endl;
}




vector<socialPlace> build_world_simple(vector<SPtype> spt,
									   vector<unsigned int> n_sp,
									   vector< probaDistrib<unsigned int> > p_size,
									   vector<individual>& indiv,
									   vector<areaUnit> auvec,
									   unsigned int seed ){
	
	/// Build a test world with individuals LINKED to social places.
	/// (use other function to make individuals PRESENT in social places)
	/// NOTE: this function is for test, it will eventually be useless...
	
	stopif( ( spt.size() != n_sp.size() ) ||
		   ( spt.size() != p_size.size() ),
		   "vectors must be same size");
	
	// number of types of social places
	unsigned long N_type_sp = spt.size();
	
	// Vector of vector of SP:
	// one vector for each type of SP
	// (all will be merged at the end)
	// 'y' : pre-sampled size
	vector< vector<socialPlace> > sp;
	vector< vector<unsigned int> > y;
	sp.resize(N_type_sp);
	y.resize(N_type_sp);
	
	ID cnt_id = 0; // id counterto make sure ID is uniques across all types and SP
	for(unsigned int t=0; t<N_type_sp; t++){
		sp[t].resize(n_sp[t]);
		y[t].resize(n_sp[t]);
		
		// Create empty (no indiv linked) social place:
		for (ID i=0; i<n_sp[t]; i++) {
			areaUnit A = auvec[rand() % auvec.size()];
			socialPlace x(A, cnt_id, spt[t]);
			sp[t][i] = x;
			cnt_id++;
		}
	}
	
	// vector that will increment the index when the social place gets full.
	// one index by type of SP. Initiated at 0 (:first index of vector).
	vector<ID> cnt_sp(N_type_sp,0);
	
	// pre-sample from the size distributions
	for(unsigned int t=0; t< N_type_sp; t++){
		unsigned int N = n_sp[t];  // total number of social places of this type
		probaDistrib<unsigned int> probD = p_size[t];  // distribution of social places' size
		// sample the sizes of each social places
		vector<unsigned int> y_t = probD.sample(N, seed+t);
		y[t] = y_t;
	}
	
	bool stoploop = false;
	vector<bool> all_sp_this_type_filled(N_type_sp, false);
	unsigned int i = 0;
	
	for(i=0; i< indiv.size() && !stoploop; i++){
		// DEBUG
		//cout << endl;
		for(ID t=0; t< N_type_sp; t++)
		{
			// DEBUG
			//			unsigned int dummy = n_sp[t];
			//			cout << "indiv_id: "<< i << "  sp_type: "<<t<< "  count:" << cnt_sp[t]
			//			<<"  y:"<<y[t][cnt_sp[t]]
			//			<< "  linked: "<< sp[t][cnt_sp[t]].n_linked_indiv() << endl;
			//
			//			if(y[t][cnt_sp[t]]==0){
			//				cout << "problem" <<endl;
			//				exit(99);
			//			}
			
			// ---------------------
			
			// Add link to current social place (which is not full)
			if(sp[t][cnt_sp[t]].n_linked_indiv() < y[t][cnt_sp[t]])
				//			   &&  cnt_sp[t] <= n_sp[t] -1)
			{
				indiv[i].set_id_sp(spt[t], sp[t][cnt_sp[t]]);
			}
			
			// If max size reached for current social place,
			// and remains social places (of same type) to be filled,
			// link this individual to the NEXT social place of the same type.
			else if(sp[t][cnt_sp[t]].n_linked_indiv() == y[t][cnt_sp[t]] &&
					cnt_sp[t] < n_sp[t]-1){
				cnt_sp[t] = cnt_sp[t] + 1;
				indiv[i].set_id_sp(spt[t], sp[t][cnt_sp[t]]);
				
				if(cnt_sp[t] % 500 == 0){
					cout << SPtype2string(spt[t])<<" : " << cnt_sp[t] <<"/"<<n_sp[t]<<endl;
				}
			}
			
			// Current social place reached its max capicity
			// and is the last one of this type.
			// Do nothing.
			else if	(sp[t][cnt_sp[t]].n_linked_indiv() == y[t][cnt_sp[t]] &&
					 cnt_sp[t] == n_sp[t]-1){
				// recorde this type is saturated:
				all_sp_this_type_filled[t] = true;
			}
			
		}
		// Test if all social places of any type reached
		// their respective size of linked indiv:
		bool tmp = 1;
		for(int t=0; t < N_type_sp; t++)
			tmp = tmp * all_sp_this_type_filled[t];
		if (tmp) {
			stoploop = true;}
	}
	stopif(i>=indiv.size(),"All individuals exhausted before filling all social places!");
	
	vector<socialPlace> spfinal = sp[0];
	
	// merge all vectors into a single one:
	for(unsigned int i=1; i<sp.size(); i++)
		spfinal.insert(spfinal.end(), sp[i].begin(), sp[i].end() );
	
	return spfinal;
}


vector<ID> at_least_one_indiv_present(const vector<socialPlace>& x){
	/// Return the IDs of the social places that have
	/// at least one indivial present
	
	vector<ID> res;
	for (ID k=0; k<x.size(); k++) {
		if (x[k].get_size()>0) res.push_back(x[k].get_id_sp());
	}
	return res;
}


void socialPlace::time_update(double dt){
	/// Update clock of all individuals in this social place
	for (unsigned int i=0; i<_indiv.size(); i++) {
		string event = _indiv[i].time_update(dt);

		if (event == "E_to_I" && _indiv[i].is_symptomatic() ) {
			_n_Is++;
			_n_E--;
		}
		if (event == "E_to_I" && !_indiv[i].is_symptomatic() ) {
			_n_Ia++;
			_n_E--;
		}
		
		if (event == "I_to_R" && _indiv[i].was_symptomatic() ) {
			_n_Is--;
			_n_R++;
		}
		if (event == "I_to_R" && !_indiv[i].was_symptomatic() ) {
			_n_Ia--;
			_n_R++;
		}
	}
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
























