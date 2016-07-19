//
//  simulation.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-06.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "simulation.h"


void Simulation::base_constructor(){
	
	_n_E = 0;
	_n_Ia = 0;
	_n_Is = 0;
	_n_R  = 0;
	_incidence = 0;
	_horizon = -999;
	
	_ts_times.clear();
	_ts_E.clear();
	_ts_Ia.clear();
	_ts_Is.clear();
	_ts_R.clear();
	
}

Simulation::Simulation(){
	base_constructor();
}


void Simulation::build_single_world(unsigned int n_indiv){
	/// Build a very simple world:
	/// one social place with 'n_indiv' individuals.
	/// USED FOR TEST
	
	cout << endl << "Building single world..."<<endl;
	
	string region1 = "Halton";
	ID id_region1 = 1;
	
	areaUnit A1(1, "singletown", id_region1, region1);
	
	vector<areaUnit> A {A1};
	
	// Schedule
	int nt = _modelParam.get_prm_uint("nt");
	vector<double> timeslice(nt, 1.0/nt); // must sum up to 1.0
	vector<SPtype> single_sed  {SP_household};
	schedule sched_worker_sed(single_sed, timeslice, "worker_sed");
	
	vector<schedule> sched {
		sched_worker_sed
	};
	
	// Disease stage durations:
	string dol_distrib = "exp";
	string doi_distrib = "exp";
	
	// individuals
	unsigned int num_indiv			= (unsigned int)(n_indiv * 3); // <-- make sure it's large enough
	vector<individual> many_indiv	= build_individuals(num_indiv,
														sched,
														dol_distrib,
														doi_distrib);
	
	// type of social places
	vector<SPtype> spt {SP_household};
	
	// populate social places with individuals (built above)
	unsigned int num_hh      = 1;

	cout << "Number of household places: " << num_hh <<endl;

	
	// W A R N I N G
	// same order as type of social places definition ('spt')
	vector<unsigned int> num_sp {num_hh};
	
	// Distribution of the size of each social place type:
	probaDistrib<unsigned int> p_hh({n_indiv},{1.0});
	
	vector<probaDistrib<unsigned int> > p_size {p_hh};
	
	// build the world:
	vector<socialPlace> W = build_world_simple(spt, num_sp, p_size, many_indiv, A);
	set_world(W);
	
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
	cout << "... test world built." << endl;
}


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
	
	// Disease stage durations:
	string dol_distrib = "exp";
	string doi_distrib = "exp";
	

	
	// individuals
	unsigned int num_indiv			= (unsigned int)(1e7 * sizereduction); // <-- make sure it's large enough
	vector<individual> many_indiv	= build_individuals(num_indiv,
														sched,
														dol_distrib,
														doi_distrib);
	
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
	set_world(W);
	
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


void Simulation::time_update(double dt){
	/// Make all relevant updates when time is advanced by 'dt'
	
	// Main timer:
	_current_time += dt;
	
	// Update individuals' clock:
	for (unsigned int k=0; k<_world.size(); k++) {
		_world[k].time_update(dt);
	}
	update_pop_count();
}


void Simulation::test(){}


void Simulation::run(){
	/// Run the simulated epidemic
    
    
    // Retrieve all model parameters:
    double p_move = _modelParam.get_prm_double("proba_move");
    bool debug_mode = _modelParam.get_prm_bool("debug_mode");
    
	if(debug_mode) cout << endl << endl << " ======= START SIMULATION ======" <<endl<<endl;
	
	_current_time = 0.0;
	
	// TO DO: CHANGE THAT, IT's UGLY AND DANGEROUS
	ID ii = at_least_one_indiv_present(_world)[0];
	vector<double> timeslice = _world[ii].get_indiv()[0].get_schedule().get_timeslice();
	// - - - - - - - - -
	
	unsigned long nts = timeslice.size();
	unsigned int k = 0;
	

	// MAIN LOOP FOR TIME
	
	for (_current_time=0.0; _current_time < _horizon; ) {
		
		if(debug_mode) check_book_keeping();
		
		unsigned int idx_timeslice = k % nts;
//		cout << "iter = " << k << " ; time slice = " << idx_timeslice;
//		cout << " ; currTime = " << _current_time <<endl;
		
		double dt = timeslice[idx_timeslice];
		
		// Actions:
		
//		cout << "BEFORE move" << endl;
//		display_split_pop_present();
		
		if(p_move>0) move_individuals_sched(idx_timeslice, p_move);
		
//		cout << "AFTER move, BEFORE transmission" << endl;
//		display_split_pop_present();
		transmission_world(dt);
//		cout << "AFTER transmissiom" << endl;
//		display_split_pop_present();

		// Update counts of population categories
		// at the 'simulation' level
		// (social places were updated during transmission):
		update_pop_count();
		
		// Record for time series:
		_ts_times.push_back(_current_time);
		_ts_incidence.push_back(_incidence);
		_ts_prevalence.push_back(_prevalence);
		_ts_S.push_back(_n_S);
		_ts_E.push_back(_n_E);
		_ts_Ia.push_back(_n_Ia);
		_ts_Is.push_back(_n_Is);
		_ts_R.push_back(_n_R);
		
		// Advance time:
		time_update(dt);
		update_pop_count();
		k++;
	}

    if(debug_mode){
        cout << endl << endl << "Simulation completed."<< endl;
        displayVector(_ts_incidence);
        cout <<"Final size: " << sumElements(_ts_incidence)<<endl;
    }
}


void Simulation::set_world(world w){
	_world = w;
	update_pop_count();
}

void Simulation::set_disease(const disease &d){
	/// Set the disease 'd' to all individuals in all social places
	for (unsigned int k=0; k<_world.size(); k++) {
		_world[k].set_disease_to_all_indiv(d);
	}
}


void Simulation::move_individuals_sched(unsigned int idx_timeslice, double proba){
	/// Move individuals across social places according to their schedule
	
	unsigned long N = _world.size();
	
	std::uniform_real_distribution<double> unif(0.0,1.0);
	
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
				double u = unif(_RANDOM_GENERATOR);
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
	
	std::uniform_real_distribution<double> unif(0.0,1.0);
	
	for (int k=0; k<_world.size(); k++)
	{
		vector<unsigned int> pos2move; // idx position to move (must be done at the end of the loop, bc vector keeps on changing size!)
		for (int i=0; i<_world[k].get_size(); i++)
		{
			// Draw the chance move will actually happen:
			double u = unif(_RANDOM_GENERATOR);
			if ( u < proba )
			{
				individual tmp = _world[k].get_indiv(i);
				
				ID id_dest = __UNDEFINED_ID;
				
				if(sptype == SP_household)	id_dest = tmp.get_id_sp_household();
				if(sptype == SP_workplace)	id_dest = tmp.get_id_sp_workplace();
				if(sptype == SP_school)	id_dest = tmp.get_id_sp_school();
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
	/// Performs transmission within the k^th social place.
	/// Returns incidence for THIS social place, during the time step 'dt'
	
	stopif(k >= _world.size(), "Asking for an inexistent social place");
	
	unsigned int inc = 0;
	
	// DEBUG
	//_world[k].displayInfo();
	// ---
	
	// Count categories of individuals:
	unsigned long N   = _world[k].get_size();
	unsigned int nIa  = _world[k].get_n_Ia();
	unsigned int nIs  = _world[k].get_n_Is();
	unsigned int nS   = _world[k].get_n_S();

	unsigned long N2 = _world[k].get_indiv().size();
	
	stopif( (nIa+nIs >N) || (nS > N) || (N != N2),
		   " DANGER : book keeping problem!");

	// Calculate number of contacts:
	
	bool homog_cont = _modelParam.get_prm_bool("homogeneous_contact");
	unsigned int nContactsDrawn;
	vector<unsigned int> pos_s;
	
	if(!homog_cont){
		double nContacts = (double)(nIa + nIs) * contact_rate * dt;
		std::poisson_distribution<> poiss(nContacts);
		
		nContactsDrawn = poiss(_RANDOM_GENERATOR);
		if (nContactsDrawn>nS) nContactsDrawn = nS;
		
		cout << "DEBUG: nContacts = " << nContacts << " ; nContactsDrawn = " <<nContactsDrawn;
		cout << " ; nS = "<< nS << endl;

		// Choose randomly the susceptibles contacted:
		// TO DO: optimize when nContactsDrawn = nS (no need to pick them, take them all!)
		pos_s = _world[k].pick_rnd_susceptibles(nContactsDrawn);

	}
	else {
		// Homogeneous contact: an infectious
		// makes contact with ALL susceptible
		// and has a probability to transmit
		// at each contact.
        double mean_n_contacts = contact_rate * dt * nS / N ;
		std::poisson_distribution<> poiss(mean_n_contacts);
//        cout << "DEBUG: mean contacts = "<< mean_n_contacts << endl;
		for (unsigned int i=0; i<(nIa+nIs); i++) {
            unsigned int inc_k =poiss(_RANDOM_GENERATOR);
			inc += inc_k;
//            cout << "DEBUG: #contact_"<<k<<" = "<< inc_k << endl;
		}
		pos_s = _world[k].pick_rnd_susceptibles(inc);
        for (unsigned int j=0; j<pos_s.size(); j++){
            _world[k].acquireDisease(pos_s[j]);
            update_pop_count();
        }
	} // end-else

	
	// DEBUG
//	_world[k].displayInfo();
//	cout << " IDs susceptible picked: ";
//	for	(int m=0; m<pos_s.size(); m++) cout<< _world[k].get_indiv()[pos_s[m]].get_id()<<", ";
//	cout<<endl;
	
	// Attempt transmission:
	// TO DO: put that in a member function of socialPlace
	
	if(!homog_cont)
	{
		std::uniform_real_distribution<double> unif(0.0,1.0);
		double p = -999, u = -999;
		
		for (unsigned int j=0; j<pos_s.size(); j++)
		{
			// Probability for THIS susceptible to acquire the disease:
			
			if (homog_cont)  p = contact_rate * dt;
			if (!homog_cont) p = _world[k].get_indiv(pos_s[j]).calc_proba_acquire_disease();
			
			u = unif(_RANDOM_GENERATOR);
			if(u < p) {
				
				//DEBUG
				//			cout << " BEFORE transmission"<<endl;
				//			_world[k].displayInfo();
				//			cout << "indiv ID_" <<_world[k].get_indiv()[pos_s[j]].get_id()<<" is going to acquire"<<endl;
				
				// Transmission!
				_world[k].acquireDisease(pos_s[j]);
				update_pop_count();
				
				//			cout << "check infected (must=1): " << _world[k].get_indiv()[pos_s[j]].is_infected() <<endl;
				//_world[k].displayInfo(); // DEBUG
				
				inc++;
				
				//DEBUG
				//			cout << "indiv ID_" <<_world[k].get_indiv()[pos_s[j]].get_id()<<" acquired!"<<endl;
				//			cout << " AFTER transmission"<<endl;
				//			_world[k].displayInfo();
			}
		}
	} // end-if-homog-cont
	return inc;
}


void Simulation::transmission_world(double timeslice){
	/// Simulates disease transmissions in the whole world (all social places)
	
	unsigned int incidence = 0;
	double cr = _modelParam.get_prm_double("contact_rate");
	
	for(unsigned int k=0; k < _world.size(); k++){
		incidence += transmission_oneSP(k, cr, timeslice);
	}
	_incidence = incidence;
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
			
			// check:
			string errmsg =  "Cannot seed in SP_ID_" + to_string(i) + " because it has less individuals than intial infections requested!";
			stopif(_world[i].get_size() < I0[cnt], errmsg);

			// seed infection in this social place:
			for(unsigned int k=0; k<I0[cnt]; k++)	_world[i].acquireDisease(k);
			_world[i].set_n_E(I0[cnt]);
			cnt++;
		}
	}
	update_pop_count();
}


void Simulation::displayInfo_indiv(){
	/// Display informations on all individuals, in all social places:
	
	for (unsigned int k=0; k<_world.size(); k++) {
		for (ID i=0; i<_world[k].get_size(); i++) {
			_world[k].get_indiv(k).displayInfo();
		}
	}
	
	
}


dcDataFrame Simulation::timeseries(){
	
	dcDataFrame df(_ts_times,"time");
	
	df.addcol("incidence", to_vector_double(_ts_incidence));
	df.addcol("prevalence", to_vector_double(_ts_prevalence));
	df.addcol("nS", to_vector_double(_ts_S));
	df.addcol("nE", to_vector_double(_ts_E));
	df.addcol("nIa", to_vector_double(_ts_Ia));
	df.addcol("nIs", to_vector_double(_ts_Is));
	df.addcol("nR", to_vector_double(_ts_R));
	
	return df;
}



void Simulation::update_pop_count(){
	/// Update the count of individuals in a all
	/// stages of the disease.
	
	unsigned long n = _world.size();
	_n_S  = 0;
	_n_E  = 0;
	_n_Ia = 0;
	_n_Is = 0;
	_n_R  = 0;
	for (ID i=0; i<n; i++) {
		_n_S  += _world[i].get_n_S();
		_n_E  += _world[i].get_n_E();
		_n_Ia += _world[i].get_n_Ia();
		_n_Is += _world[i].get_n_Is();
		_n_R  += _world[i].get_n_R();
	}
	_prevalence = _n_E + _n_Ia + _n_Is;
}




void Simulation::check_book_keeping(){
	/// Check consistency of book keeping.
	/// WARNING: FOR DEBUG ONLY, SLOWS DOWN EXECUTION!
	
	// stage S
	unsigned int nS=0;
	unsigned int nS_census=0;
	for (ID k=0; k<_world.size(); k++) {
		nS += _world[k].get_n_S();
		nS_census += _world[k].census_disease_stage("S");
	}
	bool check_S = ( (_n_S == nS) && (nS == nS_census) );
	stopif(!check_S, "Book keeping error with S stage");
	
	// stage E
	unsigned int nE=0;
	unsigned int nE_census=0;
	for (ID k=0; k<_world.size(); k++) {
		nE += _world[k].get_n_E();
		nE_census += _world[k].census_disease_stage("E");
	}
	bool check_E = ( (_n_E == nE) && (nE == nE_census) );
	
//	if(0){ // DEBUG
//		vector<ID> x = _world[0].census_disease_stage_ID("E");
//		displayVector(x);
//	}
	stopif(!check_E, "Book keeping error with E stage");
	
	
	// stage Is
	unsigned int nIs=0;
	unsigned int nIs_census=0;
	for (ID k=0; k<_world.size(); k++) {
		nIs += _world[k].get_n_Is();
		nIs_census += _world[k].census_disease_stage("Is");
	}
	bool check_Is = ( (_n_Is == nIs) && (nIs == nIs_census) );
	stopif(!check_Is, "Book keeping error with Is stage");

	// stage Ia
	unsigned int nIa=0;
	unsigned int nIa_census=0;
	for (ID k=0; k<_world.size(); k++) {
		nIa += _world[k].get_n_Ia();
		nIa_census += _world[k].census_disease_stage("Ia");
	}
	bool check_Ia = ( (_n_Ia == nIa) && (nIa == nIa_census) );
	stopif(!check_Ia, "Book keeping error with Ia stage");
	
	// stage R
	unsigned int nR=0;
	unsigned int nR_census=0;
	for (ID k=0; k<_world.size(); k++) {
		nR += _world[k].get_n_R();
		nR_census += _world[k].census_disease_stage("R");
	}
	bool check_R = ( (_n_R == nR) && (nR == nR_census) );
	stopif(!check_R, "Book keeping error with R stage");
	
}








