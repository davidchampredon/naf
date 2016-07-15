//
//  probaDistribution.h
//  naf
//
//  Created by David CHAMPREDON on 2016-07-08.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#ifndef __naf__probaDistribution__
#define __naf__probaDistribution__

#include <stdio.h>
#include <vector>
#include <random>

#include "utils.h"

using namespace std;

template <class T>
class probaDistrib {
	
	vector<T>		_value;
	vector<double>	_proba;
	
public:
	probaDistrib(){}
	probaDistrib(vector<T> value, vector<double> proba){
		stopif(value.size()!=proba.size(), "value and proba must be same size");
		_value = value;
		_proba = proba;
	}
	
	vector<T> sample(unsigned int n, unsigned int seed);
};

template <class T> vector<T> probaDistrib<T>::sample(unsigned int n, unsigned int seed){
	
	std::mt19937 gen(seed);
	std::discrete_distribution<> d(_proba.begin(), _proba.end());
	
	vector<T> s(n);
	for (int i=0; i<n; i++) s[i] = _value[d(gen)];
	return s;
}

#endif /* defined(__naf__probaDistribution__) */
