//
//  areaUnit.h
//  naf
//
//  Created by David CHAMPREDON on 2016-07-06.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#ifndef __naf__areaUnit__
#define __naf__areaUnit__


#include <stdio.h>
#include <string>
#include "dcTools.h"

using namespace std;

class areaUnit{
	

protected:
	// Features of this area unit:
	string	_name_au;
	ID		_id_au;
	double	_latitude;
	double	_longitude;
	
	// Region that includes several area units:
	ID		_id_region;
	string	_name_region;

public:
	
	// Constructors:
	areaUnit();
	areaUnit(ID id_au, string name_au, ID id_region, string regionName) {
		_id_au = id_au ;
		_name_au = name_au;
		_name_region = regionName;
		_id_region = id_region;
	}
	
	
	// Get functions:
	string	get_name_au(){return _name_au;}
	ID		get_id_au(){return _id_au;}
	double	get_latitude(){return _latitude;}
	double	get_longitude(){return _longitude;}
	ID		get_id_region(){return _id_region;}
	string	get_name_region(){return _name_region;}
	
	
	// Misceleanous:
	void displayInfo();
	void displayInfo_AU();
	
};



#endif /* defined(__naf__areaUnit__) */
