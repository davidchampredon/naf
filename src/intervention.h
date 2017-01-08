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
#include "individual.h"
#include "globalvar.h"

using namespace std;

class intervention{
    
private:
    
    float _time_start;
    float _time_end;
    
    string _name;
    string _type_intervention;
    string _type_indiv_targeted;

    float  _cvg_rate;
    float  _cvg_max_proportion;
    
    float  _efficacy;
    
    
public:
    
    intervention();
    intervention(string type_treatment,
                 string type_indiv_targeted,
                 string name,
                 float time_start, float time_end,
                 float cvg_rate, float cvg_max_prop,
                 float efficacy);

    // Set functions
    
    void set_cvg_rate(float x) {_cvg_rate = x;}
    void set_cvg_max_proportion(float x) {_cvg_max_proportion = x;}
    
    // Get functions
    
    string get_type_intervention()  const {return _type_intervention;}
    string get_type_indiv_targeted()const {return _type_indiv_targeted;}
    float get_time_start()          const {return _time_start;}
    float get_time_end()            const {return _time_end;}
    float get_cvg_rate()            const {return _cvg_rate;}
    float get_cvg_max_proportion()  const {return _cvg_max_proportion;}
    
    // Actions
    
    /** Activate intervention at the individual level according to intervention type.
     */
    void    act_on_individual(vector<individual*> x,
                              float current_time,
                              float doi_reduc_treat,
                              float imm_hum_incr,
                              float imm_cell_incr,
                              float frail_incr,
                              float vax_lag);
    
    /** Vaccinate selected symptomatic individuals. */
    void    vaccinate(vector<individual*> x,
                      float current_time,
                      float imm_hum_incr,
                      float imm_cell_incr,
                      float frail_incr,
                      float lag);
    
    /** Treat selected symptomatic individuals. */
    void    treat(vector<individual*> x, float doi_reduction);
    
    /** Instantaneously cure individual (--> doi_drawn=0). Used for debuging. */
    void    cure(vector<individual*> x);
    
    void    display_info();
    
};

#endif /* defined(__naf__intervention__) */
