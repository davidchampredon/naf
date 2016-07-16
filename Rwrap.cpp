//
//  Rwrap.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-15.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include <Rcpp.h>

using namespace Rcpp;

#include "tests.h"



// [[Rcpp::export]]
List naf_test(List params){
	
	vector<unsigned int> res;
	
	// Unpack parameters:
	unsigned int rnd_seed = params["rnd_seed"];
	cout << "DEBUG: the seed is "<<rnd_seed <<endl;
	
	// Set random seed
	_RANDOM_GENERATOR.seed(rnd_seed);
	
	// Call C++ function
	res = test_transmission();
	
	// Return R-formatted result:
	return List::create( Named("res_naf") = res);
}




// [[Rcpp::export]]
List rnd_test(List params){
	
	unsigned long N = 1e1;
	unsigned int rnd_seed = params["rnd_seed"];
	
	_RANDOM_GENERATOR.seed(rnd_seed);
	test_random();
//	std::uniform_real_distribution<double> unif(0.0,1.0);
//	cout << endl << "SEED IS: "<< rnd_seed << endl;
//	for (int i=0; i<N; i++) {
//		double y = unif(_RANDOM_GENERATOR);
//		cout << "TEST-distrib = " << y <<endl;
//	}
	cout << endl;
	
	
	
	_RANDOM_GENERATOR.seed(rnd_seed);
	test_random();
	
	return List::create( Named("res_rnd") = N);
}