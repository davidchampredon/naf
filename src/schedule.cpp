//
//  schedule.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-07.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "schedule.h"


uint pos_schedname(string name, vector<schedule> sched){
    /// Find the position of a schedule named 'name'
    /// in a vector of schedules
    
    uint pos = 0;
    
    while (sched[pos].get_name() != name ) {
        pos++;
        stopif(pos >= sched.size(),"Schedule named "+name+ " not found in vector of schedules.");
    }
    
    return pos;
}