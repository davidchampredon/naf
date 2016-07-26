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

	uint	_incidence;
	uint	_prevalence;
    
    // Number of individuals:
    uint	_n_S;   // susceptible
	uint	_n_E;   // latent stage
	uint	_n_Ia;  // infectious stage, asymptomatic
	uint	_n_Is;  // infectious stage, symptomatic
	uint	_n_R;   // recovered stage
    uint	_n_H;   // hospitalized
	
	// time series
	vector<double>	_ts_times;
	vector<uint>	_ts_incidence;
	vector<uint>	_ts_prevalence;
	vector<uint>	_ts_S;
	vector<uint>	_ts_E;
	vector<uint>	_ts_Ia;
	vector<uint>	_ts_Is;
    vector<uint>	_ts_H;
	vector<uint>	_ts_R;
    
    dcDataFrame     _ts_census_by_SP;

	
public:
	
	modelParam _modelParam;
	
	// Constructors
	void base_constructor();
	Simulation();

	// pseudo constructors:
	void build_test_world(double reduction_size);
    void build_test_2_sp(uint n_indiv);
	void build_single_world(uint n_indiv);
    void build_test_hospitalization(uint n_indiv);
	
	// Simulate
	void run();
	
	
	// Set functions
	
	void set_current_time(double t) {_current_time = t;}
	void set_world(world w);
	void set_horizon(double h) {_horizon = h;}
	void set_modelParam(modelParam mp) {_modelParam = mp;}

	
	// Get functions
	
	world	get_world() {return _world;}
	
	double			get_current_time()   const {return _current_time;}
	vector<double>	get_ts_times()       const {return _ts_times;}
	vector<uint>	get_ts_incidence()   const {return _ts_incidence;}
    dcDataFrame     get_ts_census_by_SP()const {return _ts_census_by_SP;}
	
	// Time updates
	
	void	time_update(double dt);
	void	update_pop_count();
	
	// Migration

	void move_individuals_sched(uint idx_timeslice, double proba);
	void move_individuals(const SPtype sptype, double proba);
	void move_one_individual(uint pos_indiv, ID from, ID to);
    void assign_hospital_to_individuals();
    
    // Book keeping
    
    void    check_book_keeping();
    void    define_all_id_tables();
	
	// Epidemic
	
	void	set_disease(const disease& d);
	void	seed_infection(vector<ID> id_sp, vector<uint> I0);
	uint	transmission_oneSP(uint k, double contact_rate, double dt);
	void	transmission_world(double timeslice);
	uint	prevalence();

    double          draw_contact_rate(individual* indiv, uint k);
    vector<uint>    draw_n_contacts(uint k,
                                    double dt,
                                    string infectious_type);
    
    vector< vector<uint> > draw_contacted_S(uint k,
                                            vector<uint> n_contacts,
                                            string infectious_type);
    
    double          calc_proba_transmission(individual* infectious,
                                            individual* susceptible);
    
    double          calc_proba_symptomatic(float immunity, float frailty);
    double          calc_proba_hospitalized(float frailty);
    
    vector< vector<uint> > transmission_attempts(uint k,
                                                 vector< vector<uint> > selected_S);

    uint    transmission_activation(int k,
                                    vector< vector<uint> > selected_S,
                                    vector< vector<uint> > transm_success);
    void    hospitalize_indiv(uint k, uint i);
    void    hospitalize();
    
    
	// Exports
	
	dcDataFrame		timeseries();
    
    
    // Census
    
   	uint	census_total_alive();
    void    update_ts_census_by_SP();
    
    
    // Miscelleanous

	uint	population_size();
	void	display_split_pop_present();
	void	display_split_pop_linked();
	
	void	displayInfo_indiv();
	
	void test();
	
	
	
};




#endif /* defined(__naf__simulation__) */



