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
#include "dcTools.h"
#include "probaDistribution.h"

#include "tests.h"

//#include "schedule.h"

using namespace std;


/*
 TO DO:
 - recovery from disease
 
 */



int main(int argc, const char * argv[]) {
	
    auto t0 = std::chrono::system_clock::now();
    
	double horizon = 90.0;
	
	modelParam MP;

	MP.add_prm_bool("debug_mode", true);
    
	MP.add_prm_double("dol_mean", 2.0);
	MP.add_prm_double("doi_mean", 3.0);
	MP.add_prm_double("proba_move", 0.0);
	MP.add_prm_bool("homogeneous_contact", true);
	MP.add_prm_double("contact_rate", 1);
    MP.add_prm_uint("nt", 3);
	
	unsigned int n_indiv = 1E4;
	unsigned int i0 = 2;
	
	_RANDOM_GENERATOR.seed(123);
	Simulation sim1 = test_transmission(MP,horizon,n_indiv,i0);
	
	sim1.get_world()[0].export_dcDataFrame().display();
	sim1.timeseries().display();
	
//	_RANDOM_GENERATOR.seed(123);
//	Simulation sim2 = test_transmission(MP,horizon);
//	
//	_RANDOM_GENERATOR.seed(123);
//	test_rnd_eng();
//	_RANDOM_GENERATOR.seed(123);
//	test_rnd_eng();

    // timers:
    auto t2 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = t2-t0;
    cout.precision(3);
    cout << "TOTAL TIME ELAPSED: "<< elapsed_seconds.count()/60.0 << " minutes" <<endl;
    
	return 0;
}







