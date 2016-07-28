//
//  intervention.h
//  naf
//
//  Created by David CHAMPREDON on 2016-07-27.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#ifndef __naf__intervention__
#define __naf__intervention__

#include <stdio.h>
#include <string>
#include "dcTools.h"

using namespace std;

class intervention{
    
private:
    
    float _time_start;
    float _time_end;
    
    string _name;
    string _type;

    float   _cvg_rate;
    float   _cvg_max_proportion;
    
    
public:
    
    intervention();
    intervention(string type, string name,
                 float time_start, float time_end,
                 float cvg_rate, float cvg_max_prop);

    // Set functions
    
    void set_cvg_rate(float x) {_cvg_rate = x;}
    void set_cvg_max_proportion(float x) {_cvg_max_proportion = x;}
    
    // Get functions
    
    float get_time_start()      const {return _time_start;}
    float get_time_end()        const {return _time_end;}
    float get_cvg_rate()        const {return _cvg_rate;}
    float get_cvg_max_proportion()   const {return _cvg_max_proportion;}
    
    
};

#endif /* defined(__naf__intervention__) */
