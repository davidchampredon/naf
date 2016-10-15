//
//  discrete_prob_dist.h
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

#include "dcTools.h"
#include "globalvar.h"

using namespace std;

template <class T>
class discrete_prob_dist {
    
    vector<T>		_value;
    vector<double>	_proba;
    
public:
    
    discrete_prob_dist(){}
    
    discrete_prob_dist(vector<T> value, vector<double> proba){
        _value = value;
        _proba = proba;
        check_construction();
    }
    
    
    discrete_prob_dist<T>(string file){
        dcMatrix A = distributionFromFile(file);
        vector<double> val0  = A.extractColumn(0);
        vector<double> proba = A.extractColumn(1);
        vector<T> val(val0.begin(), val0.end());
        _value = val;
        _proba = proba;
        check_construction();
    }
    
    unsigned long size() const { return _value.size();}
    vector<T> get_value() const {return _value;}
    
    vector<T> sample(uint n, uint seed);
    
    void display(){
        vector<double> val(_value.begin(), _value.end());
        dcMatrix A (val);
        A.addColVector(_proba);
        A.display();
    }
    
    void check_construction(){
        stopif(_value.size()!=_proba.size(), "Value and proba must be same size.");
        double s = sumElements(_proba);
        stopif( abs(s-1.0) > 1E-6 , "Probabilities do not sum up to 1.0. They currently sum up to " + to_string(s));
    }
};


template <class T> vector<T> discrete_prob_dist<T>::sample(uint n, uint seed){
    /// Draw 'n' samples from this distribution
    
    std::discrete_distribution<> d(_proba.begin(), _proba.end());
    
    vector<T> s(n);
    for (int i=0; i<n; i++) s[i] = _value[d(_RANDOM_GENERATOR)];
    return s;
}





#endif /* defined(__naf__probaDistribution__) */
