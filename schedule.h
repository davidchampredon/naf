//
//  schedule.h
//  naf
//
//  Created by David CHAMPREDON on 2016-07-07.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#ifndef __naf__schedule__
#define __naf__schedule__

#include <stdio.h>
#include <vector>
#include <numeric>
#include <cmath>

//#include "socialPlace.h"
#include "dcTools.h"

using namespace std;

class schedule{
	
	vector<double>	_timeslice;
	vector<SPtype>	_sp_type;
	string			_name;
	
public:
	
	schedule(){}
	
	schedule(vector<double> timeslice){
		double sum = accumulate(timeslice.begin(), timeslice.end(), 0.0);
		stopif(fabs(sum - 1.0)>1e-6, "Vector 'timeslice' must sum up to 1.0 (currently sum is " + to_string(sum) + ")");
		_timeslice = timeslice;
	}

	schedule(vector<SPtype> sp_type, vector<double> timeslice, string name){
		double sum = accumulate(timeslice.begin(), timeslice.end(), 0.0);
		stopif(fabs(sum - 1.0)>1e-6, "Vector 'timeslice' must sum up to 1.0 (currently sum is " + to_string(sum) + ")");

		_sp_type	= sp_type;
		_timeslice	= timeslice;
		_name		= name;
	}
	
	// Set functions:
	
	void set_sp_type(vector<SPtype> sp_type) {_sp_type = sp_type;}
	
	
	// Get functions:
	vector<SPtype>	get_sp_type() {return _sp_type;}
	string			get_name() {return _name;}
	vector<double>	get_timeslice(){return _timeslice;}
};

uint pos_schedname(string name, vector<schedule> sched);

#endif /* defined(__naf__schedule__) */
