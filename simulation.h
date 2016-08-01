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
#include "intervention.h"

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
    uint    _n_D;   // dead
    
    uint    _n_treated;
    uint    _n_vaccinated;
	
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
    vector<uint>    _ts_D;
    vector<uint>    _ts_n_treated;
    vector<uint>    _ts_n_vaccinated;

    dcDataFrame     _ts_census_by_SP;

    vector<intervention> _intervention;
    
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
	uint	transmission_oneSP(uint k, double dt);
	void	transmission_world(double timeslice);
	uint	prevalence();

    double          draw_contact_rate(individual* indiv, uint k);
    vector<uint>    draw_n_contacts(uint k,
                                    double dt,
                                    string infectious_type);
    
    vector< vector<uint> > draw_contacted_S(uint k,
                                            vector<uint> n_contacts,
                                            string infectious_type);
    
    double calc_proba_transmission(individual* infectious,
                                            individual* susceptible);
    
    double calc_proba_symptomatic(float immunity, float frailty);
    double calc_proba_hospitalized(float frailty);
    double calc_proba_death(float frailty);
    
    uint    transmission_process(uint k, double dt, string infectious_type);
    
    vector< vector<uint> > transmission_attempts(uint k,
                                                 vector< vector<uint> > selected_S,
                                                 string infectious_type);

    void    transmission_wiw(int k,
                             vector< vector<uint> > selected_S,
                             vector< vector<uint> > transm_success,
                             string infectious_type);
    
    uint    transmission_activation(int k,
                                    vector< vector<uint> > selected_S,
                                    vector< vector<uint> > transm_success);
    void    hospitalize_indiv(uint k, uint i);
    void    hospitalize();
    void    discharge_hospital(uint idx_timeslice);
    void    death_hospital();

    // Interventions
    
    void    add_intervention(intervention x) {_intervention.push_back(x);}
    void    activate_interventions(ID id_sp, double dt);
    void    update_immunity_frailty();
    bool    at_least_one_vaccination_intervention();
    vector<individual*> draw_targeted_individuals(uint i_intervention,
                                                  ID id_sp,
                                                  double dt);
    
    
    // Census
    
   	uint	census_total_alive();
    void    update_ts_census_by_SP();

    
	// Exports
	
	dcDataFrame		timeseries();
    
    
    
    // Miscelleanous

	uint	population_size();
	void	display_split_pop_present();
	void	display_split_pop_linked();
	
	void	displayInfo_indiv();
	
	void test();
	
	
	
};




#endif /* defined(__naf__simulation__) */



