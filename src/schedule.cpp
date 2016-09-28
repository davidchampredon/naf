//
//  schedule.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-07.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "schedule.h"


double schedule::name_to_double(){
    /// Translate the name (string) of the schedule
    /// into a double. Value of double returned
    /// has no meaning, just a numerical way
    /// to differentiate schedule type. Helper function.
    
    bool found = false;
    double res = 0;
    
    if      (_name=="worker_sed")   res = 1;
    else if (_name=="student")      res = 2;
    else if (_name=="unemployed")   res = 3;
    
    found = (res>0);
    stopif(!found, "Double translation for schedule name <"+_name+"> unknonwn.");
    
    return res;
}


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