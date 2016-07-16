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
	
	res = test_transmission();
	
	return List::create( Named("res_naf") = res);
	
}