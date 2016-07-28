//
//  intervention.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-27.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "intervention.h"



intervention::intervention(){
    _time_start = __UNDEFINED_FLOAT;
    _time_end = __UNDEFINED_FLOAT;
    
    _cvg_rate = 0.0;
    _cvg_max_proportion = 0.0;
    _name = "UNDEFINED INTERVENTION";
    _type = "UNDEFINED INTERVENTION";
    
}

intervention::intervention(string type, string name,
                           float time_start, float time_end,
                           float cvg_rate, float cvg_max_prop){
    _type = type;
    _name = name;
    _time_start = time_start;
    _time_end = time_end;
    _cvg_rate = cvg_rate;
    _cvg_max_proportion = cvg_max_prop;
}