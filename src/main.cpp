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
#include <chrono>

#include <random>

#include "individual.h"
#include "socialPlace.h"
#include "simulator.h"
#include "dcTools.h"
#include "discrete_prob_dist.h"

#include "tests.h"


using namespace std;


int main(int argc, const char * argv[]) {
    
    auto t0 = std::chrono::system_clock::now();
    
    system("pwd");
    
    try{
        main_run_test();
    }
    catch (...) {
        std::cout << endl <<  "   ~~~~~ Standard C++ exception caught ~~~~~ "<<endl <<endl;
    }

    // Previous tests:
    // main_test_naf(); 
    // main_test_SEIR_vs_ODE();
    // main_test_move_2_sp();
    
    // timers:
    auto t2 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = t2-t0;
    cout.precision(3);
    cout << "TOTAL TIME ELAPSED: "<< elapsed_seconds.count()/60.0 << " minutes" <<endl;
    
    return 0;
}







