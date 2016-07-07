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
//	double cr = _modelParam.get_prm_double("contact_rate");
	//move_individuals(SP_household, p);
	//transmission_oneSP(1, cr, 1.0);
	move_individuals_sched(2, p);
}


void Simulation::run(){
	/// Run the simulated epidemic
	
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
		cout << "k = " << k << " ; time slice = " << idx_timeslice << " ; currTime = " << _current_time <<endl;
		
		// Actions:
		move_individuals_sched(idx_timeslice, p);
		transmission_world(timeslice[idx_timeslice]);

		// Record for time series:
		_ts_times.push_back(_current_time);
		_ts_incidence.push_back(_current_incidence);
		
		// Advance time:
		_current_time += timeslice[idx_timeslice];
		k++;
	}
	cout << endl << endl << " DEBUG::simulation completed."<< endl;
}



void Simulation::move_individuals_sched(unsigned int idx_timeslice, double proba){
	
	/// Move individuals across social places according to their schedule
	
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

				// Retrieve the social place this individual is supposed
				// to go for this timeslice, according to its schedule
				SPtype sptype = tmp.get_schedule().get_sp_type()[idx_timeslice];
				
//				cout << " DEBUG: sptype: " << SPtype2string(sptype) <<endl;
				
				// Retrieve the actual destination:
				
				ID id_dest = __UNDEFINED_ID;
				if(sptype == SP_household)	id_dest = tmp.get_id_sp_household();
				if(sptype == SP_workplace)	id_dest = tmp.get_id_sp_workplace();
				if(sptype == SP_school)		id_dest = tmp.get_id_sp_school();
				if(sptype == SP_other)		id_dest = tmp.get_id_sp_other();
				if(sptype == SP_hospital)	id_dest = tmp.get_id_sp_hospital();
				if(sptype == SP_pubTransp)	id_dest = tmp.get_id_sp_pubTransp();
				
//				cout << " DEBUG: SP ID: " << id_dest <<endl;
				
				// if destination not defined or
				// if indiv is already here, do nothing
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
	
	stopif(k >= _world.size(), "asking for an inexistent social place");
	
	// number of infected
	// TO DO: distinguish INFECTED vs INFECTIOUS
	unsigned int n = _world[k].get_prevalence();
	// number of susceptible in this social place:
	unsigned int allpop = (unsigned int)(_world[k].get_size());
	stopif(allpop < n, "BOOK KEEPING PROBLEM!");
	unsigned int nS = allpop - n;
	
	unsigned int inc = 0;
	
	if (n > 0 && nS > 0){
		
		// Calculate total number of contacts:
		unsigned int nContacts = (unsigned int)(contact_rate * n * dt);
		
		cout << " DEBUG TRANSMISSION: nContacts: " << nContacts << endl;
		
		
		
		// Susceptible candidates for transmission:
		vector<unsigned int> susc = _world[k].pick_rnd_susceptibles(nContacts);
		
		
		cout << " DEBUG TRANSMISSION: #susc: " << susc.size() << endl;
		
		// Calculate transmission based on susceptible features:
		
		for (unsigned int i=0; i<susc.size(); i++) {
			
			// Transmission Probability
			
			double imm = _world[k].get_indiv()[susc[i]].get_immunity();
			double fra = _world[k].get_indiv()[susc[i]].get_frailty();
			double pt = (1-imm) * fra ;
			
			cout << " DEBUG TRANSMISSION: pt: " << pt << endl;
			
			
			// Draw event:
			double u = (double)(rand())/RAND_MAX;
			if( u < pt) {
				_world[k].get_indiv()[susc[i]].acquireDisease();
				_world[k].increase_prevalence();
				inc ++;
			}
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