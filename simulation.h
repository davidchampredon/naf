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
#include "globalvar.h"
#include "dcDataFrame.h"

using world = vector<socialPlace>;


class Simulation{

protected:
	
	world	_world;
	double	_horizon;
	double	_current_time;
	
	unsigned int	_incidence;
	unsigned int	_prevalence;
	unsigned int	_n_S;   // susceptible
	unsigned int	_n_E;   // latent stage
	unsigned int	_n_Ia;  // infectious stage, asymptomatic
	unsigned int	_n_Is;  // infectious stage, symptomatic
	unsigned int	_n_R;   // recovered stage
	
	// time series
	vector<double>			_ts_times;
	vector<unsigned int>	_ts_incidence;
	vector<unsigned int>	_ts_prevalence;
	vector<unsigned int>	_ts_S;
	vector<unsigned int>	_ts_E;
	vector<unsigned int>	_ts_Ia;
	vector<unsigned int>	_ts_Is;
	vector<unsigned int>	_ts_R;

	
public:
	
	modelParam _modelParam;
	
	// Constructors
	void base_constructor();
	Simulation();

	// pseudo constructors:
	void build_test_world(double reduction_size);
	void build_single_world(unsigned int n_indiv);
	
	// Simulate
	void run();
	
	
	// Set functions
	
	void set_current_time(double t) {_current_time = t;}
	void set_world(world w);
	void set_horizon(double h) {_horizon = h;}
	void set_modelParam(modelParam mp) {_modelParam = mp;}

	
	// Get functions
	
	world	get_world() {return _world;}
	
	double					get_current_time()	{return _current_time;}
	vector<double>			get_ts_times()		{return _ts_times;}
	vector<unsigned int>	get_ts_incidence()	{return _ts_incidence;}
	
	// Time updates
	
	void	time_update(double dt);
	void	update_pop_count();
	
	// Migration

	void move_individuals_sched(unsigned int idx_timeslice, double proba);
	void move_individuals(const SPtype sptype, double proba);
	void move_one_individual(unsigned int k,unsigned int i, const SPtype sptype);
	
    
    // Book keeping
    
    void			check_book_keeping();
    void            define_all_id_tables();
	
	// Epidemic
	
	void			set_disease(const disease& d);
	void			seed_infection(vector<ID> id_sp, vector<unsigned int> I0);
	unsigned int	transmission_oneSP(unsigned int k, double contact_rate, double dt);
	void			transmission_world(double timeslice);
	unsigned int	prevalence();
	

	// Exports
	
	dcDataFrame		timeseries();
    
    
	// Miscelleanous

	unsigned int	census_total_alive();
	
	

	
	unsigned int	population_size();
	void			display_split_pop_present();
	void			display_split_pop_linked();
	
	void			displayInfo_indiv();
	
	void test();
	
	
	
};


#endif /* defined(__naf__simulation__) */



