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
	double	_current_time;
	
	unsigned int	_current_incidence;
	
	
	// time series
	vector<double>			_ts_times;
	vector<unsigned int>	_ts_incidence;
	
	
	
public:
	
	modelParam _modelParam;
	
	// Constructors
	
	Simulation();
	Simulation(world w, double h) {_world = w; _horizon = h;}
	
	// Simulate
	void run();
	
	
	// Set functions
	
	void set_current_time(double t) {_current_time = t;}
	
	// Get functions
	
	world	get_world() {return _world;}
	
	double					get_current_time()	{return _current_time;}
	vector<double>			get_ts_times()		{return _ts_times;}
	vector<unsigned int>	get_ts_incidence()	{return _ts_incidence;}
	
	
	// Migration

	void move_individuals_sched(unsigned int idx_timeslice, double proba);
	void move_individuals(const SPtype sptype, double proba);
	void move_one_individual(unsigned int k,unsigned int i, const SPtype sptype);
	
	
	// Epidemic
	
	void seed_infection(vector<ID> id_sp, vector<unsigned int> I0);
	unsigned int transmission_oneSP(unsigned int k, double contact_rate, double dt);
	void transmission_world(double timeslice);
	unsigned int prevalence();
	
	// Miscelleanous

	unsigned int census_total_alive();
	unsigned int population_size();
	void displayPopulationSplit();
	void test();
	
	
	
};

// DELETE??
//inline void acquireDisease(individual& x){
//	x.acquireDisease();
//}


#endif /* defined(__naf__simulation__) */



