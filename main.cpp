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
    
    
     main_test_naf();
    
    // Previous tests:
    // main_test_SEIR_vs_ODE();
    // main_test_move_2_sp();
    
    // timers:
    auto t2 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = t2-t0;
    cout.precision(3);
    cout << "TOTAL TIME ELAPSED: "<< elapsed_seconds.count()/60.0 << " minutes" <<endl;
    
    return 0;
}







