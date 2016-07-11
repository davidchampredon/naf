//
//  simulation.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-06.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "simulation.h"




void Simulation::build_test_world(double sizereduction){
	
	/// Build a test world with individuals LINKED and PRESENT to social places.
	/// NOTE: this function is for test, it will eventually be useless...
	
	cout << endl << "Building test world..."<<endl;
	
	string region1 = "Halton";
	ID id_region1 = 1;
	
	areaUnit A1(1, "Oakville", id_region1, region1);
	areaUnit A2(2, "Burlington", id_region1, region1);
	areaUnit A3(3, "Hamilton", id_region1, region1);
	
	vector<areaUnit> A {A1,A2,A3};
	
	// Schedule
	vector<double> timeslice {2.0/24, 8.0/24, 2.0/24, 12.0/24}; // must sum up to 1.0
	vector<SPtype> worker_sed  {SP_pubTransp, SP_workplace, SP_pubTransp, SP_household};
	vector<SPtype> student     {SP_pubTransp, SP_school,    SP_pubTransp, SP_household};
	
	schedule sched_worker_sed(worker_sed, timeslice, "worker_sed");
	schedule sched_student(student, timeslice, "student");
	
	//	vector<SPtype> worker_trav {SP_pubTransp, SP_workplace, SP_other, SP_household};
	//	schedule sched_worker_trav(worker_trav, timeslice, "worker_trav");
	
	vector<schedule> sched {
		sched_worker_sed,
		sched_student
	};
	
	// individuals
	unsigned int num_indiv			= (unsigned int)(1e7 * sizereduction); // <-- make sure it's large enough
	vector<individual> many_indiv	= build_individuals(num_indiv, sched);
	
	// type of social places
	vector<SPtype> spt {
		SP_pubTransp,
		SP_workplace,
		SP_household,
		SP_school
	};
	
	// populate social places with individuals (built above)
	unsigned int num_pubTr   = (unsigned int)(2000 * sizereduction);
	unsigned int num_biz     = (unsigned int)(180e3 * sizereduction);
	unsigned int num_hh      = (unsigned int)(1500e3 * sizereduction);
	unsigned int num_school  = (unsigned int)(2600 * sizereduction);
	
	cout << "Number of public transportations places: " << num_pubTr <<endl;
	cout << "Number of business places: " << num_biz <<endl;
	cout << "Number of household places: " << num_hh <<endl;
	cout << "Number of school places: " << num_school <<endl;
	
	// W A R N I N G
	// same order as type of social places definition ('spt')
	vector<unsigned int> num_sp {
		num_pubTr,
		num_biz,
		num_hh,
		num_school
	};
	
	// Distribution of the size of each social place type:
	
	probaDistrib<unsigned int> p_pubTransp({20,30,60},{0.6,0.3,0.1});
	probaDistrib<unsigned int> p_workPlace({3,7,15,30,75,200},{0.6,0.15,0.15,0.07,0.02,0.01});
	probaDistrib<unsigned int> p_hh({1,2,3,4,5,6,7,8},{0.23, 0.34, 0.16, 0.15, 0.06, 0.03, 0.02, 0.01});
	probaDistrib<unsigned int> p_school({250,500,750,1000,1250,1500},{0.60,0.25,0.10,0.03, 0.01,0.01});
	
	vector<probaDistrib<unsigned int> > p_size {
		p_pubTransp,
		p_workPlace,
		p_hh,
		p_school
	};
	
	// build the world:
	vector<socialPlace> W = build_world_simple(spt, num_sp, p_size, many_indiv, A);
	
	// assign to class members:
	_world = W;
	
	// initial population: Move everyone to its household!
	unsigned long N = _world.size();
	for (int k=0; k<N; k++)
	{
		if(_world[k].get_type()==SP_household)
		{
			unsigned int n_linked_k = _world[k].n_linked_indiv();
			for (unsigned int i=0; i<n_linked_k; i++)
			{
				ID curr_indiv_ID = _world[k].get_linked_indiv_id()[i];
				individual tmp = get_indiv_with_ID(curr_indiv_ID, many_indiv);
				
				bool debugcode = false;
				if(debugcode){
					// Checks (remove for better speed)
					ID id_hh = tmp.get_id_sp_household();
					stopif(id_hh == __UNDEFINED_ID, "at least one individual has no linked household!");
					stopif(id_hh != k, "Not consistent linkage!");
					_world[id_hh].add_indiv(tmp);
					// -----
				}
				
				// Faster version:
				if(!debugcode) _world[k].add_indiv(tmp);
			}
		}
	}
	cout << "... test world built."<<endl;
}





void Simulation::test(){}


void Simulation::run(){
	/// Run the simulated epidemic
	
	cout << endl << endl << " ======= START SIMULATION ======" <<endl<<endl;
	
	_current_time = 0.0;
	
	// TO DO: CHANGE THAT, IT's UGLY AND DANGEROUS
	ID ii = at_least_one_indiv_present(_world)[0];
	vector<double> timeslice = _world[ii].get_indiv()[0].get_schedule().get_timeslice();
	// - - - - - - - - -
	unsigned long nts = timeslice.size();
	unsigned int k = 0;
	
	// Retrieve all model parameters:
	double p = _modelParam.get_prm_double("proba_move");
	
	
	// MAIN LOOP FOR TIME
	
	for (_current_time=0.0; _current_time < _horizon; ) {
		
		unsigned int idx_timeslice = k % nts;
		cout << "iter = " << k << " ; time slice = " << idx_timeslice << " ; currTime = " << _current_time <<endl;
		
		// Actions:
		
//		cout << "BEFORE move" << endl;
//		display_split_pop_present();
		
		move_individuals_sched(idx_timeslice, p);
		
//		cout << "AFTER move, BEFORE transmission" << endl;
//		display_split_pop_present();

		transmission_world(timeslice[idx_timeslice]);
		
//		cout << "AFTER transmissiom" << endl;
//		display_split_pop_present();

		// Record for time series:
		_ts_times.push_back(_current_time);
		_ts_incidence.push_back(_current_incidence);
		
		// Advance time:
		_current_time += timeslice[idx_timeslice];
		k++;
		
//		cout << "DEBUG: total prev: " << prevalence() << endl;
//		cout << "DEBUG: census alive: " << census_total_alive() << endl;
//		cout << "DEBUG: population_size: " << population_size()  << endl;
	}
	cout << endl << endl << "Simulation completed."<< endl;
	displayVector(_ts_incidence);
	cout <<"Final size: " << sumElements(_ts_incidence)<<endl;
}



void Simulation::move_individuals_sched(unsigned int idx_timeslice, double proba){
	
	/// Move individuals across social places according to their schedule
	
	
	unsigned long N = _world.size();
	
	for (int k=0; k<N; k++)
	{
		for (unsigned int i=0; i<_world[k].get_size(); i++)
		{
			// Retrieve its actual destination
			ID id_dest = _world[k].find_dest(i, idx_timeslice);
			
			if ( (id_dest!=k) && (id_dest!=__UNDEFINED_ID) )
			{
				// take a copy of the individual
				individual tmp = _world[k].get_indiv(i);
				
				// Draw the chance move will actually happen:
				// TO DO: make this proba individual-dependent
				double u = (double)(rand()) / RAND_MAX;
				if ( u < proba ){
					
					//DEBUG
//					cout << endl << "BEFORE move" <<endl;
//					_world[k].displayInfo();
					// ------
					
					// add individual at destination
					_world[id_dest].add_indiv(tmp);
					// remove this individual (in i^th position in '_indiv' vector) from here
					_world[k].remove_indiv(i);
					
					
					//DEBUG
//					cout << endl << "AFTER move" <<endl;
//					_world[k].displayInfo();
					// ------
					

				}
			}
		}
	}
	
}


void Simulation::move_individuals(const SPtype sptype, double proba){
	
	/// Move individuals across social places
	
	for (int k=0; k<_world.size(); k++)
	{
		vector<unsigned int> pos2move; // idx position to move (must be done at the end of the loop, bc vector keeps on changing size!)
		for (int i=0; i<_world[k].get_size(); i++)
		{
			// Draw the chance move will actually happen:
			double u = (double)(rand()) / RAND_MAX;
			if ( u < proba )
			{
				individual tmp = _world[k].get_indiv(i);
				
				ID id_dest = __UNDEFINED_ID;
				
				if(sptype == SP_household)	id_dest = tmp.get_id_sp_household();
				if(sptype == SP_workplace)	id_dest = tmp.get_id_sp_workplace();
				if(sptype == SP_school)		id_dest = tmp.get_id_sp_school();
				if(sptype == SP_other)		id_dest = tmp.get_id_sp_other();
				if(sptype == SP_hospital)	id_dest = tmp.get_id_sp_hospital();
				if(sptype == SP_pubTransp)	id_dest = tmp.get_id_sp_pubTransp();
				
				if(id_dest != __UNDEFINED_ID &&
				   id_dest != k)
				{
					// add indiv to destination
					_world[id_dest].add_indiv(tmp);
					// record its position for future deletion
					pos2move.push_back(i);
				}
			}
		}
		// remove from this SP the individuals that moved:
		if(pos2move.size()>0) _world[k].remove_indiv(pos2move);
	}
}


unsigned int Simulation::transmission_oneSP(unsigned int k,
											double contact_rate,
											double dt){
	/// Performs transmission within the k^th social place. Returns incidence
	
	stopif(k >= _world.size(), "Asking for an inexistent social place");
	
	unsigned int inc = 0;
	
	// Count categories of individuals:
	
	// DEBUG
	//_world[k].displayInfo();
	// ---
	
	unsigned long N = _world[k].get_size();
	unsigned int nI = _world[k].get_prevalence();
	stopif(nI > N, " DANGER : book keeping problem!");
	unsigned int nS = (unsigned int)(N - nI);
	
	// Calculate number of contacts:
	unsigned int nContacts = (unsigned int)(nI * contact_rate * dt);
	if (nContacts>nS) nContacts = nS;
	
	// Choose randomly the susceptibles contacted:
	// TO DO: optimize when nContact = nS (no need to pick them, take them all!)
	vector<unsigned int> pos_s = _world[k].pick_rnd_susceptibles(nContacts);
	
	
	// DEBUG
//	_world[k].displayInfo();
//	cout << " IDs susceptible picked: ";
//	for	(int m=0; m<pos_s.size(); m++) cout<< _world[k].get_indiv()[pos_s[m]].get_id()<<", ";
//	cout<<endl;
	
	// Attempt transmission:
	// TO DO: put that in a member function of socialPlace
	for (unsigned int j=0; j<pos_s.size(); j++)
	{
		double p = _world[k].get_indiv(pos_s[j]).calc_proba_acquire_disease();
		double u = (double)(rand())/RAND_MAX;
		if(u < p) {
			
			//DEBUG
//			cout << " BEFORE transmission"<<endl;
//			_world[k].displayInfo();
//			cout << "indiv ID_" <<_world[k].get_indiv()[pos_s[j]].get_id()<<" is going to acquire"<<endl;
			
			// Transmission!
			_world[k].acquireDisease(pos_s[j]);
			
//			cout << "check infected (must=1): " << _world[k].get_indiv()[pos_s[j]].is_infected() <<endl;
			//_world[k].displayInfo(); // DEBUG
			
			inc++;
			
			//DEBUG
//			cout << "indiv ID_" <<_world[k].get_indiv()[pos_s[j]].get_id()<<" acquired!"<<endl;
//			cout << " AFTER transmission"<<endl;
//			_world[k].displayInfo();
		}
	}
	return inc;
}


void Simulation::transmission_world(double timeslice){
	/// Simulates disease transmissions in the whole world (all social places)
	
	unsigned int incidence = 0;
	double cr = _modelParam.get_prm_double("contact_rate");
	
	for(unsigned int k=0; k < _world.size(); k++){
		incidence += transmission_oneSP(k, cr, timeslice);
	}
	_current_incidence = incidence;
}


unsigned int Simulation::census_total_alive(){
	/// Counts all individuals that are alive
	unsigned int cnt = 0;
	for(int k=0; k<_world.size(); k++) cnt += _world[k].census_alive();
	return cnt;
}


unsigned int Simulation::prevalence(){
	/// Prevalence in the whole world
	
	unsigned int cnt = 0;
	for (int k=0; k<_world.size(); k++) cnt += _world[k].get_prevalence();
	return cnt;
}


unsigned int Simulation::population_size(){
	
	unsigned int s = 0;
	for(int i=0; i<_world.size(); i++) s+=_world[i].get_size();
	return s;
}


void Simulation::display_split_pop_present(){
	
	cout<<endl<<"------------"<<endl;
	unsigned int s = 0;
	unsigned int p = 0;
	for(int i=0;i<_world.size();i++){
		cout << "sp_"<<i<<" : present = " << _world[i].get_size()<<" (prev="<<_world[i].get_prevalence()<<")";
		cout << "\t ["<< SPtype2string(_world[i].get_type())<<"]" <<endl;
		s += _world[i].get_size();
		p += _world[i].get_prevalence();
	}
	cout<< "Total population = "<< s << " (prev="<<p<<")"<<endl;
	cout<<"------------"<<endl;
	
	
}

void Simulation::display_split_pop_linked(){
	cout<<endl<<"------------"<<endl;
	unsigned int s = 0;
	for(int i=0;i<_world.size();i++){
		cout << "sp_"<<i<<" : Linked indiv = " << _world[i].n_linked_indiv();
		cout << "\t ["<< SPtype2string(_world[i].get_type())<<"]" <<endl;
		s += _world[i].n_linked_indiv();
	}
	cout<< "Total linked = "<< s << endl;
	cout<<"------------"<<endl;
}

void Simulation::seed_infection(vector<ID> id_sp, vector<unsigned int> I0){
	/// Seed infection in specified socialplaces, with specified initial number of infectious indiv
	
	stopif(id_sp.size() != I0.size(), "vectors must be same size");
	
	ID cnt = 0;
	for(ID i=0; i<_world.size(); i++){
		if (_world[i].get_id_sp() == id_sp[cnt]) {
			stopif(_world[i].get_size() < I0[cnt], "Cannot seed in SP_ID_" + to_string(i) + " because it has less individuals than intial infections requested!");
			for(unsigned int k=0; k<I0[cnt]; k++)	_world[i].acquireDisease(k);
			cnt++;
		}
	}
}


