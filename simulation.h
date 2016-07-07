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
#include "modelParam.h"


using world = vector<socialPlace>;


class Simulation{

protected:
	
	world	_world;
	double	_horizon;
	double	_timestep;
	
	
public:
	
	modelParam _modelParam;
	
	// Constructors
	
	Simulation();
	Simulation(world w, double h) {_world = w; _horizon = h;}
	
	// Get functions
	
	world get_world() {return _world;}
	
	
	// Migration

	void move_individuals_sched(unsigned int idx_timeslice, double proba);
	void move_individuals(const SPtype sptype, double proba);
	void move_one_individual(unsigned int k,unsigned int i, const SPtype sptype);
	
	
	
	// Epidemic
	

	unsigned int transmission_oneSP(unsigned int k, double contact_rate, double dt);
	
	unsigned int prevalence();
	
	// Miscelleanous

	unsigned int census_total_alive();
	void test();
	
	
	
};

// DELETE??
inline void acquireDisease(individual& x){
	x.acquireDisease();
}


#endif /* defined(__naf__simulation__) */



