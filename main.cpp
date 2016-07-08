//
//  main.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-04.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <limits.h>

#include <random>

#include "individual.h"
#include "socialPlace.h"
#include "simulation.h"
#include "utils.h"

//#include "schedule.h"

using namespace std;


/*
 TO DO:
 - create large number of individuals, with age distribution
 - create large number of social places, associate individuals
 -
 
 */



int main(int argc, const char * argv[]) {
	
	
//
//	vector<double> y {0.2, 0.5, 0.3};
//	vector<SPtype> w {SP_school, SP_household, SP_hospital};
//	schedule sched(y);
//	sched.set_sp_type(w);
//
//	exit(99);
	
	string region1 = "Halton";
	ID id_region1 = 1;
	
	areaUnit A1(1, "Oakville", id_region1, region1);
	areaUnit A2(2, "Burlington", id_region1, region1);
	areaUnit A3(3, "Hamilton", id_region1, region1);
	
	A1.displayInfo();
	
	vector<areaUnit> A {A1,A2,A3};
	
	// Schedule
	vector<double> timeslice {0.3, 0.5, 0.2};
	vector<SPtype> worker_sed {SP_pubTransp, SP_workplace, SP_household};
	vector<SPtype> worker_trav {SP_pubTransp, SP_household, SP_school};
	
	schedule sched_worker_sed(worker_sed, timeslice, "worker_sed");
	schedule sched_worker_trav(worker_trav, timeslice, "worker_trav");
	vector<schedule> sched {sched_worker_sed,sched_worker_trav};
	
	
	double horiz = 20;
	unsigned int N_sp		= 100;
	unsigned int N_indiv	= 5000;
	
	world W = build_random(N_sp, A);
	populate_random_with_indiv(W, N_indiv, sched);
	displayPopulationSize(W);
	

	Simulation sim(W, horiz);
	
	sim.seed_infection({0,1}, {2,1});
	
	sim.displayPopulationSplit();
	
	sim._modelParam.add_prm_double("proba_move", 0.90);
	sim._modelParam.add_prm_double("contact_rate", 3.0);

	sim.run();
	
	return 0;
	/*
	// Define social places:
	socialPlace sp1(A1,0, SP_school);
	socialPlace sp2(A2,1, SP_household);
	socialPlace sp3(A3,2, SP_workplace);
	socialPlace sp4(A3,3, SP_pubTransp);
	
//	sp1.displayInfo();
//	sp2.displayInfo();
//	sp3.displayInfo();
//	sp4.displayInfo();
	
	vector<socialPlace> spvec;
	spvec.push_back(sp1);
	spvec.push_back(sp2);
	spvec.push_back(sp3);
	spvec.push_back(sp4);
	
	
	
	
	// Individuals
	vector<individual> indivvec;
	ID n_indiv = 24;
	
	for(int i=0; i<n_indiv; i++){
		double age = rand() % 90 + 1;
		
		individual tmp(i, age);
		tmp.set_id_sp_school(sp1);
		tmp.set_id_sp_household(sp2);
		tmp.set_id_sp_workplace(sp3);
		tmp.set_id_sp_pubTransp(sp4);
		
		tmp.set_immunity(0.0);
		tmp.set_frailty(1.0);
		
		tmp.set_schedule(sched_worker_sed);
		if(i>10) tmp.set_schedule(sched_worker_trav);
		
		indivvec.push_back(tmp);
		//tmp.displayInfo();
	}
	
	indivvec[5].acquireDisease();
	indivvec[8].acquireDisease();
	indivvec[2].forget_id_sp_household();
	
	// assign individuals to SP
	for (int i=0; i<indivvec.size(); i++) {
		int sp_idx = rand() % spvec.size();
		spvec[sp_idx].add_indiv(indivvec[i]);
	}
	
	
	
	double horizon = 20;
	Simulation sim(spvec, horizon);
	
	sim._modelParam.add_prm_double("proba_move", 0.90);
	sim._modelParam.add_prm_double("contact_rate", 3.0);
	
//	cout <<"total prev0 = "<< sim.prevalence()<<endl;
	
//	sim.displayPopulationSplit();
//	
//	sim.test();
//	
//	sim.displayPopulationSplit();
	
	sim.get_world()[0].displayInfo();
	sim.get_world()[1].displayInfo();
	sim.get_world()[2].displayInfo();
	sim.get_world()[3].displayInfo();
	
//	cout <<"total pop = "<< sim.census_total_alive()<<endl;
//	cout <<"total prev = "<< sim.prevalence()<<endl;
	
//	cout << "BEFORE RUN" <<endl;
//	sim.displayPopulationSplit();
	
	sim.displayPopulationSplit();
	sim.run();
	sim.displayPopulationSplit();
	cout << "prevalence= "<<sim.prevalence()<<endl;
//	cout << "AFTER RUN" <<endl;
//	sim.displayPopulationSplit();

	
	
	//	// Locations
	//	vector<string> locname;
	//	locname.push_back("Burlington");
	//	locname.push_back("Hamilton");
	//	locname.push_back("Oakville");
	//
	//	vector<location> loc;
	//	for(int i=0; i<locname.size(); i++){
	//		location tmp(i, locname[i]);
	//		loc.push_back(tmp);
	//	}
	//
	//	// Individuals
	//	vector<individual> indiv;
	//	ID n_indiv = 24;
	//	for(int i=0; i<n_indiv; i++){
	//		double age = rand() % 90 + 1;
	//		individual tmp(i, age);
	//		indiv.push_back(tmp);
	//		tmp.displayInfo();
	//	}
	//
	//	// Households
	//	vector<household> hhvec;
	//
	//	vector<int> hhloc = {0,1,2,1};
	//	for(int i=0; i<hhloc.size(); i++){
	//		household tmp(i,loc[hhloc[i]]);
	//		hhvec.push_back(tmp);
	//	}
	//
	//	for (int i=0; i<hhvec.size(); i++){
	//		hhvec[i].displayInfo();
	//	}
	//
	//	vector<ID> idvec;
	//
	//	for(int i=0;i<=11;i++) idvec.push_back(i);
	//	hhvec[0].populate_household(indiv, idvec);
	//
	//	idvec.clear();
	//	for(int i=12;i<=16;i++) idvec.push_back(i);
	//	hhvec[1].populate_household(indiv, idvec);
	//	idvec.clear();
	//	for(int i=17;i<=19;i++) idvec.push_back(i);
	//	hhvec[2].populate_household(indiv, idvec);
	//	idvec.clear();
	//	for(int i=20;i<=23;i++) idvec.push_back(i);
	//	hhvec[3].populate_household(indiv, idvec);
	//
	//
	//	hhvec[3].displayInfo();
	//	indiv[21].displayInfo();
	//
	//
	//	vector<double> age;
	//	vector<double> proba;
	//	age.push_back(2); proba.push_back(0.2);
	//	age.push_back(11); proba.push_back(0.2);
	//	age.push_back(22); proba.push_back(0.2);
	//	age.push_back(33); proba.push_back(0.2);
	//	age.push_back(44); proba.push_back(0.2);
	//	age.push_back(55); proba.push_back(0.2);
	//
	//
	//	vector<double> proba_loc;
	//	proba_loc.push_back(0.6);
	//	proba_loc.push_back(0.3);
	//	proba_loc.push_back(0.1);
	//
	//	ID popsize = 30;
	//	population P;
	//	P.create_indiv(popsize);
	//	P.assign_location(loc, proba_loc, 1234);
	//	P.assign_age(age, proba, 1234);
	//	P.displayInfo();
	//
	//
	*/
	
	return 0;
}
