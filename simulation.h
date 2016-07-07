//
//  simulation.h
//  naf
//
//  Created by David CHAMPREDON on 2016-07-06.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#ifndef __naf__simulation__
#define __naf__simulation__

#include <stdio.h>

#include "individual.h"
#include "socialPlace.h"


using world = vector<socialPlace>;


// DELETE??
inline void acquireDisease(individual& x){
	x.acquireDisease();
}


void move_indiv(const SPtype sptype, world& spvec, double proba);



unsigned int transmission(socialPlace& sp, double contact_rate, double dt);


#endif /* defined(__naf__simulation__) */



