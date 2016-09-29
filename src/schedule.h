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
	vector<SPtype>	get_sp_type()           const {return _sp_type;}
    SPtype          get_sp_type(uint i)     const
    {
        stopif(i >= _sp_type.size(), "Schedule "+_name+": _sp_type["+to_string(i)+"] not defined, _sp_type.size()="+to_string(_sp_type.size()));
        return _sp_type[i];
    }
    
	string			get_name()              const {return _name;}
	vector<double>	get_timeslice()         const {return _timeslice;}
    
    // Miscellaneous
    /**
     *  Translate the name (string) of the schedule
    *  into a double. Value of double returned
    *  has no meaning, just a numerical way
    *  to differentiate schedule type. Helper function.
     */
    double          name_to_double();
    
};


/**
 * Convert a string to a social place type object.
 */
SPtype to_SPtype(string x);
/**
 * Convert a vector string to a social place type object.
 */
vector<SPtype> to_SPtype(vector<string> x);

/**
 * Find the position of a schedule named 'name'
 * in a vector of schedules.
 */
uint pos_schedname(string name, vector<schedule> sched);


/**
 * Build a vector of schedules from a 'string' description (typically read from file).
 * Format is:
 * 
 * sched_name_1    |  sched_name_2  | ...
 *
 * ---------------------------------------
 *
 * "SP_household"  |  "SP_work"     | ...
 *
 * "SP_other"      |  "SP_pubTransp"| ...
 *
 *  ...            |   ...          | ...
 */
vector<schedule> build_all_schedules(const vector<vector<string> > &x,
                                     vector<string> sched_name,
                                     vector<double> timeslice);






#endif /* defined(__naf__schedule__) */
