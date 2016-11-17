//
//  shapeFct.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-08-26.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "shapeFct.h"
#include <math.h>


double linear_interpol(double x, double x1, double y1, double x2, double y2){
    return (x-x1)/(x2-x1)*y2 + (x-x2)/(x1-x2)*y1;
}


float frailty_mean(float age,
                   float f0,
                   float fmin,
                   float agemin,
                   float agepivot,
                   float fpivot,
                   float powerChild)
{
    float res = -999.99;
    float age_infant = 2.0;
    
    if (age <= age_infant){
        res = (f0-fmin)/pow(age_infant,powerChild) * pow(age_infant-age,powerChild) + fmin;
    }
    else if (age_infant < age && age <= agemin){
        res = fmin;
    }
    else if (agemin < age){
        float b = 1/fpivot - 1;
        float alpha =  -log(1/fmin/b-1/b)/(agemin-agepivot);
        res = 1/(1+b*exp(-alpha*(age-agepivot)));
    }
    return res;
}



double immunity_humoral(float age,
                        float agezero,
                        float baseline,
                        float p)
{
    return baseline * ( pow(pow(agezero,p)-pow(age,p),1/p) )/agezero;
}



double immunity_cellular(float age,
                         float imm_max,
                         float slope,
                         float pivot){
    
    return imm_max/(1+exp(-slope*(age/pivot-1)));
}





