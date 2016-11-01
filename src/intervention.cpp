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
    _time_end   = __UNDEFINED_FLOAT;
    
    _cvg_rate = 0.0;
    _cvg_max_proportion = 0.0;
    _name = "UNDEFINED INTERVENTION";
    _type_intervention   = "UNDEFINED INTERVENTION";
    _type_indiv_targeted = "UNDEFINED INTERVENTION";
}

intervention::intervention(string type_treatment,
                           string type_indiv_targeted,
                           string name,
                           float time_start, float time_end,
                           float cvg_rate, float cvg_max_prop){
    _type_intervention = type_treatment;
    _type_indiv_targeted = type_indiv_targeted;
    _name = name;
    _time_start = time_start;
    _time_end = time_end;
    _cvg_rate = cvg_rate;
    _cvg_max_proportion = cvg_max_prop;
}


void intervention::treat(vector<individual*> x,
                                 float doi_reduction){
    /// Treat selected symptomatic individuals.
    
    for (uint i = 0; i<x.size(); i++) {
        x[i]->receive_treatment(doi_reduction);
    }
}

void intervention::vaccinate(vector<individual*> x,
                             float current_time,
                             float imm_hum_incr,
                             float imm_cell_incr,
                             float frail_incr,
                             float vaxlag){
    /// Vaccinate selected symptomatic individuals.
    
    for (uint i = 0; i<x.size(); i++) {
        x[i]->receive_vaccine(current_time, imm_hum_incr,imm_cell_incr, frail_incr, vaxlag);
    }
}



void intervention::cure(vector<individual *> x){
    /// Instantaneously cure individual (--> doi_drawn=0)
    /// Used for debuging.
    
    for (uint i = 0; i<x.size(); i++) {
        x[i]->receive_cure();
    }
}


void intervention::act_on_individual(vector<individual*> x,
                                     float current_time,
                                     float doi_reduc_treat,
                                     float imm_hum_incr,
                                     float imm_cell_incr,
                                     float frail_incr,
                                     float vax_lag){
    /// Activate intervention at the individual level
    /// according to intervention type.

    bool found = false;
    
    if (_type_intervention == "treatment") {
        treat(x,doi_reduc_treat);
        found=true;
    }
    
    else if (_type_intervention == "vaccination") {
        vaccinate(x,current_time,imm_hum_incr,imm_cell_incr,frail_incr,vax_lag);
        found = true;
    }
    
    else if (_type_intervention == "cure") {
        cure(x);
        found=true;
    }
    
    stopif(!found,"Intervention type unknown: " + _type_intervention);
}



void intervention::display_info(){
    
    cout << endl << "Intervention info:"<<endl;
    tabcout(" Name", _name);
    tabcout(" Type", _type_intervention);
    tabcout(" Targeted individuals", _type_indiv_targeted);
    tabcout(" Start intervention", _time_start);
    tabcout(" End intervention", _time_end);
    tabcout(" Coverage rate", _cvg_rate);
    tabcout(" Coverage max proportion", _cvg_max_proportion);
    cout << endl;
}



