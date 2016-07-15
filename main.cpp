//
//  main.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-04.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <limits.h>
#include <ctime>

#include <random>

#include "individual.h"
#include "socialPlace.h"
#include "simulation.h"
#include "utils.h"
#include "probaDistribution.h"

#include "tests.h"

//#include "schedule.h"

using namespace std;


/*
 TO DO:
 - recovery from disease
 
 */



int main(int argc, const char * argv[]) {
	
	test_transmission();
	
	return 0;
}
