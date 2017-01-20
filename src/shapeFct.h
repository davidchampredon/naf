//
//  shapeFct.h
//  naf
//
//  Created by David CHAMPREDON on 2016-08-26.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#ifndef __naf__shapeFct__
#define __naf__shapeFct__

#include <stdio.h>

#endif /* defined(__naf__shapeFct__) */


double linear_interpol(double x, double x1, double y1, double x2, double y2);


/** Defines mean frailty
 */
float frailty_mean(float age,
                   float f0,
                   float agepivot,
                   float slope1,
                   float slope2);

/**
 * Defines humoral immunity index for an individual of age 'age'
 */
double immunity_humoral(float age, float agezero, float baseline, float p);


/**
 * Defines cellular immunity index for an individual of age 'age'
 */
double immunity_cellular(float age, float imm_max, float slope, float pivot);





