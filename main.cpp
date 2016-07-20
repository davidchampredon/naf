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
    
//    vector<int> x1 {5,3,0,7,0,2,3,9,3};
//    vector<int> x2 {-4,1,2,8};
//    vector<int> x3 {0,0,0,1,0};
//    vector< vector<int> > z;
//    z.push_back(x1);
//    z.push_back(x2);
//    z.push_back(x3);
//    vector<int> yy = melt(z);
//    displayVector(yy);
    
    
    double horizon = 90.0;
    
    modelParam MP;
    
    MP.add_prm_bool("debug_mode", true);
    
    MP.add_prm_double("dol_mean", 2.0);
    MP.add_prm_double("doi_mean", 3.0);
    MP.add_prm_double("proba_move", 0.0);
    MP.add_prm_bool("homogeneous_contact", true);
    MP.add_prm_double("contact_rate", 1.0);
    MP.add_prm_uint("nt", 3);
    
    unsigned int n_indiv = 1e3;
    unsigned int i0 = 5;
    
    _RANDOM_GENERATOR.seed(123);
    Simulation sim1 = test_transmission(MP,horizon,n_indiv,i0);
    
    sim1.get_world()[0].export_dcDataFrame().display();
    sim1.timeseries().display();
    
    // timers:
    auto t2 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = t2-t0;
    cout.precision(3);
    cout << "TOTAL TIME ELAPSED: "<< elapsed_seconds.count()/60.0 << " minutes" <<endl;
    
    return 0;
}







