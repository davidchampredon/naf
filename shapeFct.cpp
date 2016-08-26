//
//  shapeFct.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-08-26.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "shapeFct.h"
#include <math.h>


float frailty_mean(float age,
                   float f0,
                   float fmin,
                   float agemin,
                   float agepivot,
                   float fpivot,
                   float powerChild){
    
    /// Defines mean frailty
    
    float res = -999.99;
    
    if (age<=agemin){
        res = (f0-fmin)/pow(agemin,powerChild) * pow(agemin-age,powerChild) + fmin;
    }
    else{
        float b = 1/fpivot - 1;
        float alpha =  -log(1/fmin/b-1/b)/(agemin-agepivot);
        res = 1/(1+b*exp(-alpha*(age-agepivot)));
    }
    return res;
}



