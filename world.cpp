//
//  world.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-04.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "world.h"



void world::displayInfo(){
	for (int i=0; i<_location.size(); i++)
		_location[i].displayInfo();
}