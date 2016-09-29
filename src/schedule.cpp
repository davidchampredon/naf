//
//  schedule.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-07.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "schedule.h"


double schedule::name_to_double(){
    
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
    uint pos = 0;
    while (sched[pos].get_name() != name ) {
        pos++;
        stopif(pos >= sched.size(),"Schedule named "+name+ " not found in vector of schedules.");
    }
    return pos;
}



SPtype to_SPtype(string x)
{
    SPtype res = SP_MAX;
    bool found = false;
    if      (x=="SP_household") {res = SP_household; found = true;}
    else if (x=="SP_workplace") {res = SP_workplace; found = true;}
    else if (x=="SP_school")    {res = SP_school;    found = true;}
    else if (x=="SP_other")     {res = SP_other;     found = true;}
    else if (x=="SP_hospital")  {res = SP_hospital;  found = true;}
    else if (x=="SP_pubTransp") {res = SP_pubTransp; found = true;}
    
    stopif(!found,"Cannot convert string <" +x+ "> to SPtype.");
    return res;
}


vector<SPtype> to_SPtype(vector<string> x)
{
    stopif(x.size()==0, "Empty vector.");
    vector<SPtype> res;
    for (uint i=0; i<x.size(); i++) {
        res.push_back(to_SPtype(x[i]));
    }
    return res;
}



vector<schedule> build_all_schedules(const vector<vector<string> > &sptype,
                                     vector<string> sched_name,
                                     vector<double> timeslice)
{
    size_t n = sptype.size();
    vector<schedule> schedvec(n);
    
    stopif(n != sched_name.size(),"Number of schedule names does not match schedules description.");
    
    for (size_t i=0; i<n; i++) {
        stopif(sptype[i].size() != timeslice.size(),
               "Schedule description and time slice not consistent.");
        
        vector<SPtype> x = to_SPtype(sptype[i]);
        
        schedule tmp_sched(x, timeslice, sched_name[i]);
        schedvec.push_back(tmp_sched);
    }
    return schedvec;
}



























