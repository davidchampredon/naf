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
#include "household.h"
#include "population.h"
#include "world.h"


using namespace std;


/*
 
 */



int main(int argc, const char * argv[]) {

	
	// Locations
	vector<string> locname;
	locname.push_back("Burlington");
	locname.push_back("Hamilton");
	locname.push_back("Oakville");
	
	vector<location> loc;
	for(int i=0; i<locname.size(); i++){
		location tmp(i, locname[i]);
		loc.push_back(tmp);
	}
	
	// Individuals
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
	
	vector<double> age;
	vector<double> proba;
	age.push_back(2); proba.push_back(0.2);
	age.push_back(11); proba.push_back(0.2);
	age.push_back(22); proba.push_back(0.2);
	age.push_back(33); proba.push_back(0.2);
	age.push_back(44); proba.push_back(0.2);
	age.push_back(55); proba.push_back(0.2);
	
	
	vector<double> proba_loc;
	proba_loc.push_back(0.6);
	proba_loc.push_back(0.3);
	proba_loc.push_back(0.1);
	
	ID popsize = 30;
	population P;
	P.create_indiv(popsize);
	P.assign_location(loc, proba_loc, 1234);
	P.assign_age(age, proba, 1234);
	P.displayInfo();
	
	
//	world W(loc);
//	W.displayInfo();

	
	
	return 0;
}
