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


float frailty_mean(float age,
                   float f0,
                   float fmin,
                   float agemin,
                   float agepivot,
                   float fpivot,
                   float powerChild);