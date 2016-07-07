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



using namespace std;


/*
 
 */



int main(int argc, const char * argv[]) {
	
	
//	vector<int> x {5,8,9,12};
//	
//	auto pos = distance(x.begin(), find(x.begin(),x.end(),8));
//	cout<<"TEST="<<pos<<endl;
//	
//	x.erase(x.begin()+ pos);
//	
//	exit(99);
	
	string region1 = "Halton";
	ID id_region1 = 1;
	
	areaUnit A1(1, "Oakville", id_region1, region1);
	areaUnit A2(2, "Burlington", id_region1, region1);
	areaUnit A3(3, "Hamilton", id_region1, region1);
	
	A1.displayInfo();
	
	
	// Define social places:
	socialPlace sp1(A1,0, SP_school);
	socialPlace sp2(A2,1, SP_household);
	socialPlace sp3(A3,2, SP_workplace);
	
	sp1.displayInfo();
	sp2.displayInfo();
	sp3.displayInfo();
	
	vector<socialPlace> spvec;
	spvec.push_back(sp1);
	spvec.push_back(sp2);
	spvec.push_back(sp3);
	
	
	// Individuals
	vector<individual> indivvec;
	ID n_indiv = 24;
	
	for(int i=0; i<n_indiv; i++){
		double age = rand() % 90 + 1;
		
		individual tmp(i, age);
		tmp.set_id_sp_school(sp1);
		tmp.set_id_sp_household(sp2);
		tmp.set_id_sp_workplace(sp3);
		
		tmp.set_immunity(0.0);
		tmp.set_frailty(1.0);
		
		indivvec.push_back(tmp);
		tmp.displayInfo();
	}
	
	indivvec[5].acquireDisease();
	indivvec[8].acquireDisease();
	indivvec[2].forget_id_sp_household();
	
	// assign individuals to SP
	for (int i=0; i<indivvec.size(); i++) {
		int sp_idx = rand() % 3;
		spvec[sp_idx].add_indiv(indivvec[i]);
	}
	
	spvec[0].displayInfo();
	spvec[1].displayInfo();
	spvec[2].displayInfo();
	
	cout << " - - - MOVE - - - "<<endl;
	
	move_indiv(SP_workplace, spvec, 0.90);
	
	spvec[0].displayInfo();
	spvec[1].displayInfo();
	spvec[2].displayInfo();
	
	unsigned int inc = transmission(spvec[2], 2.5 , 1.0);
	
	cout << " INCIDENCE: " << inc << endl;
	
	spvec[2].displayInfo();
	
	exit(99);
	
	
	vector<individual> sel_indiv;
	sel_indiv.push_back(indivvec[2]);
	sel_indiv.push_back(indivvec[4]);
	sel_indiv.push_back(indivvec[6]);
	
	sp1.add_indiv(sel_indiv);
	
	indivvec[5].acquireDisease();
	sp1.add_indiv(indivvec[5]);
	sp1.displayInfo();
//	indiv[4].displayInfo();
//	sp1.get_indiv()[1].displayInfo();
	indivvec[5].displayInfo();
	
	cout << "REMOVE" <<endl;
	sp1.remove_indiv(indivvec[4]);
	sp1.displayInfo();
	
	sp1.get_indiv()[0].acquireDisease();
	sp1.displayInfo();
	
	acquireDisease(sp1.get_indiv()[0]);
	sp1.displayInfo();
	
	
	sp2.add_indiv(indivvec[20]);
	sp2.add_indiv(indivvec[21]);
	sp2.displayInfo();

	cout << " TEST VECTOR"<<endl;
	
	
	
	spvec[0].displayInfo();
	spvec[1].displayInfo();
	
	
	
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
	
	
	return 0;
}
