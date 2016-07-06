//
//  world.h
//  naf
//
//  Created by David CHAMPREDON on 2016-07-04.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#ifndef __naf__world__
#define __naf__world__

#include <stdio.h>
#include <vector>

#include "location.h"

class world{
public:
	vector<location> _location;
	
	world(vector<location> vect_loc) {_location = vect_loc;}
	
	void displayInfo();
	
};

#endif /* defined(__naf__world__) */
