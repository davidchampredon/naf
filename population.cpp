//
//  population.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-04.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include <random>
#include "population.h"



//
//void population::create_indiv(ID totalNumberOfIndiv){
//	/// Create initial individuals
//	
//	_indiv.clear();
//	for(int i=0; i<totalNumberOfIndiv; i++){
//		individual tmp(i,0.0);
//		_indiv.push_back(tmp);
//	}
//	_size = totalNumberOfIndiv;
//}
//
//
//void population::assign_location(vector<location> loc, vector<double> proba, unsigned int seed){
//	/// Assign location to all individuals based on discrete probability distribution
//	
//	stopif(loc.size() != proba.size(), "number of locations and probabilities not same size!");
//	stopif(abs(sumElements(proba) - 1.0) > 1e-6, "Probabilities do not sum to 1.0!");
//	
//	vector<ID> locsize;
//	for (int i=0; i< loc.size(); i++) {
//		locsize.push_back( (unsigned int)(proba[i]*_size) );
//	}
//	ID c=0;
//	for(ID L=0; L<locsize.size(); L++){
//		if (L==0) {
//			for(ID j=0; j<locsize[0]; j++){
//				_indiv[j]._current_location = loc[L];
//				c++;
//			}
//		}
//		if(L>0){
//			ID c2 = c;
//			for (ID j=c2; j<c2+locsize[L]; j++) {
//				_indiv[j]._current_location = loc[L];
//				c++;
//			}
//			c2 = c;
//		}
//	}
//	
//}
//
//
//
//void population::assign_age(vector<double> age, vector<double> proba, unsigned int seed){
//	/// Assign ages to all individuals based on discrete probability distribution
//	
//	stopif(age.size() != proba.size(), "Ages and probabilities not same size!");
//	std::mt19937 gen(seed);
//	std::discrete_distribution<> d(proba.begin(),proba.end());
//	
//	for(int i=0; i<_size; i++) {
//		_indiv[i]._age = age[d(gen)];
//	}
//}
//
//
//void population::displayInfo(){
//	
//	for (ID i=0; i<_indiv.size(); i++) {
//		_indiv[i].displayInfo();
//	}
//	
//}