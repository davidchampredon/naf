//
//  simulation.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-06.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "simulation.h"


void Simulation::test(){
	double p = _modelParam.get_prm_double("proba_move");
	//move_individuals(SP_household, p);
	
	cout << endl <<  " - - - BEFORE MOVE - - - "<<endl;
	displayPopulationSplit();
	
	move_individuals_sched(1, p);
	
	cout << endl <<  " - - - AFTER MOVE - - - "<<endl;
	displayPopulationSplit();

	double cr = _modelParam.get_prm_double("contact_rate");
	transmission_oneSP(2, cr, 2.0);
	
	cout << endl <<  " - - - AFTER TRANSMISSION - - - "<<endl;
	displayPopulationSplit();
}


void Simulation::run(){
	/// Run the simulated epidemic
	
	cout << endl << endl << " ======= START SIMULATION ======" <<endl<<endl;
	
	_current_time = 0.0;
	
	// CHANGE THAT, IT's UGLY AND DANGEROUS
	vector<double> timeslice = _world[0].get_indiv()[0].get_schedule().get_timeslice();
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
//		displayPopulationSplit();
		
		move_individuals_sched(idx_timeslice, p);
		
//		cout << "AFTER move, BEFORE transmission" << endl;
//		displayPopulationSplit();

		transmission_world(timeslice[idx_timeslice]);
		
//		cout << "AFTER transmissiom" << endl;
//		displayPopulationSplit();

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
	cout << endl << endl << " DEBUG::simulation completed."<< endl;
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
				individual tmp = _world[k].get_indiv()[i];
				
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

//
//void Simulation::move_individuals_sched(unsigned int idx_timeslice, double proba){
//	
//	/// Move individuals across social places according to their schedule
//	
//	for (int k=0; k<_world.size(); k++)
//	{
//		
////		cout << "DEBUG: sp #"<<k<< " popsize_before_move = "<< _world[k].get_size()<<endl;
//		
//		vector<unsigned int> pos2move; // idx position to move (must be done at the end of the loop, bc vector keeps on changing size!)
//		for (int i=0; i<_world[k].get_size(); i++)
//		{
//			// Draw the chance move will actually happen:
//			double u = (double)(rand()) / RAND_MAX;
//			if ( u < proba )
//			{
//				individual tmp = _world[k].get_indiv()[i];
//
//				// Retrieve the social place this individual is supposed
//				// to go for this timeslice, according to its schedule
//				SPtype sptype = tmp.get_schedule().get_sp_type()[idx_timeslice];
//				
////				cout << " DEBUG: sptype: " << SPtype2string(sptype) <<endl;
//				
//				// Retrieve the actual destination:
//				
//				ID id_dest = __UNDEFINED_ID;
//				if(sptype == SP_household)	id_dest = tmp.get_id_sp_household();
//				if(sptype == SP_workplace)	id_dest = tmp.get_id_sp_workplace();
//				if(sptype == SP_school)		id_dest = tmp.get_id_sp_school();
//				if(sptype == SP_other)		id_dest = tmp.get_id_sp_other();
//				if(sptype == SP_hospital)	id_dest = tmp.get_id_sp_hospital();
//				if(sptype == SP_pubTransp)	id_dest = tmp.get_id_sp_pubTransp();
//				
////				cout << " DEBUG: SP ID: " << id_dest <<endl;
//				
//				// if destination not defined or
//				// if indiv is already here, do nothing
//				if(id_dest != __UNDEFINED_ID &&
//				   id_dest != k)
//				{
////					cout <<" DEBUG: indiv#"<<tmp.get_id() << " move from sp#" << k << " to sp#" << id_dest <<endl;
//					
////					cout <<"  DEBUG: popsize before sp#" << k << " = " << _world[k].get_size() <<endl;
////					cout <<"  DEBUG: popsize before sp#" << id_dest << " = " << _world[id_dest].get_size() <<endl;
//					// add indiv to destination
//					_world[id_dest].add_indiv(tmp);
//					// record its position for future deletion
//					pos2move.push_back(i);
//					
////					cout <<"  DEBUG: popsize after sp#" << k << " = " << _world[k].get_size() <<endl;
////					cout <<"  DEBUG: popsize after sp#" << id_dest << " = " << _world[id_dest].get_size() <<endl;
//					
//				}
//			}
//		}
//		// remove from this SP the individuals that moved:
//		if(pos2move.size()>0) {
//			cout << " DEBUG: sp_"<< k <<endl;
//			_world[k].remove_indiv(pos2move);
//		}
//		
////		cout << "DEBUG: sp #"<<k<< " popsize_after_move = "<< _world[k].get_size()<<endl;
//	}
//}


// DELETE?
void Simulation::move_one_individual(unsigned int k,unsigned int i, const SPtype sptype){
	/// Move only ONE individual, the i^th in the k^th social place, to its associated social place type
	
	individual tmp = _world[k].get_indiv()[i];
	
	ID id_dest = __UNDEFINED_ID;
	
	if(sptype == SP_household)	id_dest = tmp.get_id_sp_household();
	if(sptype == SP_workplace)	id_dest = tmp.get_id_sp_workplace();
	if(sptype == SP_school)		id_dest = tmp.get_id_sp_school();
	if(sptype == SP_other)		id_dest = tmp.get_id_sp_other();
	if(sptype == SP_hospital)	id_dest = tmp.get_id_sp_hospital();
	if(sptype == SP_pubTransp)	id_dest = tmp.get_id_sp_pubTransp();
	
	if(id_dest != __UNDEFINED_ID &&
	   id_dest != k){
		// add indiv to destination
		_world[id_dest].add_indiv(tmp);
		// and remove it from this SP
		_world[k].remove_indiv(i);
	}
}
// ---- DELETE ?


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
				individual tmp = _world[k].get_indiv()[i];
				
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
		double p = _world[k].get_indiv()[pos_s[j]].calc_proba_acquire_disease();
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


//unsigned int Simulation::transmission_oneSP(unsigned int k,
//											double contact_rate,
//											double dt){
//	/// Performs transmission within the k^th social place. Returns incidence
//	
//	stopif(k >= _world.size(), "asking for an inexistent social place");
//
//	unsigned int inc = 0;
//	
//	// number of infected
//	// TO DO: distinguish INFECTED vs INFECTIOUS
//	//(ie, currently SIR, should be at least SEIR)
//	unsigned int nI = _world[k].get_prevalence();
//	
//	// number of susceptible in this social place:
//	unsigned int allpop = (unsigned int)(_world[k].get_size());
//	if (allpop==0) return 0;
//	
//	stopif(allpop < nI, "BOOK KEEPING PROBLEM!");
//	unsigned int nS = allpop - nI;
//	
//	cout << " - - - - " << endl;
//	cout << k << " DEBUG TRANSMISSION: allpop: " << allpop << endl;
//	cout << k << " DEBUG TRANSMISSION: nInf: " << nI ;
//	displayVector(_world[k].id_infected_bruteforce());
//	cout << endl;
//	cout << k << " DEBUG TRANSMISSION: nS: " << nS << endl;
//	
//	
//	if (nI > 0 && nS > 0){
//		
//		// Calculate total number of contacts:
//		unsigned int nContacts = (unsigned int)(contact_rate * nI * dt);
//		if (nContacts > nS) nContacts = nS;
//		
//		cout << k << " DEBUG TRANSMISSION: nContacts: " << nContacts << endl;
//		
//		// Susceptible candidates for transmission:
//		vector<unsigned int> susc ;
//		if (nContacts < nS) susc = _world[k].pick_rnd_susceptibles(nContacts);
//		if (nContacts == nS) {for (int j=0;j<nS;j++) susc.push_back(j);}
//		
//		cout << k << " DEBUG TRANSMISSION: #susc: " << susc.size();
//		cout << " ; IDs: ";
//		for(int i=0; i<susc.size();i++) cout<<_world[k].get_indiv()[susc[i]].get_id()<<"; ";
//		cout << endl;
//		
//		// Calculate transmission based on susceptible features:
//		
//		for (unsigned int i=0; i<susc.size(); i++) {
//			
//			// Transmission Probability
//			
//			double imm = _world[k].get_indiv()[susc[i]].get_immunity();
//			double fra = _world[k].get_indiv()[susc[i]].get_frailty();
//			double pt = (1-imm) * fra ;
//			
//			//cout << " DEBUG TRANSMISSION: pt: " << pt << endl;
//			
//			
//			// Draw event:
//			double u = (double)(rand())/RAND_MAX;
//			if( u < pt) {
//				_world[k].get_indiv()[susc[i]].acquireDisease();
//				_world[k].increase_prevalence();
//				inc ++;
//			}
//		}
//	}
//	return inc;
//}


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



void Simulation::displayPopulationSplit(){
	
	cout<<endl<<"------------"<<endl;
	unsigned int s = 0;
	unsigned int p = 0;
	for(int i=0;i<_world.size();i++){
		cout << "sp_"<<i<<" : " << _world[i].get_size()<<" (prev="<<_world[i].get_prevalence()<<")"<<endl;
		s += _world[i].get_size();
		p += _world[i].get_prevalence();
	}
	cout<< "total = "<< s << " (prev="<<p<<")"<<endl;
	cout<<"------------"<<endl;
	
	
}

void Simulation::seed_infection(vector<ID> id_sp, vector<unsigned int> I0){
	/// Seed infection in specified socialplaces, with specified initial number of infectious indiv
	
	stopif(id_sp.size() != I0.size(), "vectors must be same size");
	
	ID cnt = 0;
	for(ID i=0; i<_world.size(); i++){
		if (_world[i].get_id_sp() == id_sp[cnt]) {
			stopif(_world[i].get_size()==0, "Cannot seed in SP_ID_" + to_string(i) + " because it is empty!");
			for(unsigned int k=0; k<I0[cnt]; k++)	_world[i].acquireDisease(k);
			cnt++;
		}
	}
}