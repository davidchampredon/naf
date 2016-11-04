//
//  simulator.h
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
#include "build_world.h"

using world = vector<socialPlace>;


class Simulator{

protected:
	
	world	_world;
	double	_horizon;
	double	_current_time;
    double  _start_time;

    uint    _initial_prevalence;
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
	
    vector<socialPlace*>  _sp_other;  // keep track of 'other' social places
    
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

    vector<double>  _ts_census_sp_time;
    vector<uint>    _ts_census_sp_id;
    vector<string>  _ts_census_sp_type;
    vector<uint>    _ts_census_sp_nS;
    vector<uint>    _ts_census_sp_nE;

    vector<intervention> _intervention;

    /**
     * Calculate the contact rate ratio based on features of individual and social place.
     */
    double  select_contact_rate_ratio(double age, SPtype sp_type);
    
    // trackers (mostly used for debugging)
    vector<uint>    _track_n_contacts;
    vector<uint>    _track_n_contacts_time;
    vector<uint>    _track_n_contacts_uid;
    
    vector< vector<float> > _wiw_ages;  // Nx2 matrix where col[1]=age infector, col[2]=age infectee
    
    dcMatrix        _contactAssort;
    
    /**
     * Set the 'i0' initial infectious individuals
     * in randomly chosen social places.
     */
    void initial_infections(uint i0);
    
public:
	
	modelParam _modelParam;
	
	// Constructors
	void base_constructor();
	Simulator();

	// Pseudo constructors:
    
    /**
     * Create STOCHASTICALLY a world from probability distributions inputs.
     */
    void create_world(vector<areaUnit> AU,
                      vector<discrete_prob_dist<uint> > D_size_hh,     // Households sizes
                      vector< vector<discrete_prob_dist<uint> > > pr_age_hh,  // Age distribution inside households
                      
                      vector<discrete_prob_dist<uint> > D_size_wrk,
                      vector<discrete_prob_dist<uint> > D_size_pubt,
                      vector<discrete_prob_dist<uint> > D_size_school,
                      vector<discrete_prob_dist<uint> > D_size_hosp,
                      vector<discrete_prob_dist<uint> > D_size_other,
                      vector<uint> n_hh,
                      vector<uint> n_wrk,
                      vector<uint> n_pubt,
                      vector<uint> n_school,
                      vector<uint> n_hosp,
                      vector<uint> n_other,
                      float unemployed_prop,
                      vector<schedule> sched);
    
    /**
     * Create a world from deterministic vector inputs.
     */
    void create_world_det(vector<areaUnit> AU,
                          vector<vector<uint> > size_hh,     // Households sizes
                          vector< vector<discrete_prob_dist<uint> > > pr_age_hh,  // Age distribution inside households
                          
                          vector<vector<uint> > size_wrk,
                          vector<vector<uint> > size_pubt,
                          vector<vector<uint> > size_school,
                          vector<vector<uint> > size_hosp,
                          vector<vector<uint> > size_other,

                          float unemployed_prop,
                          vector<schedule> sched);
    
    
	void build_test_world(double reduction_size);
    void build_test_2_sp(uint n_indiv);
	void build_single_world(uint n_indiv);
    void build_test_hospitalization(uint n_indiv);
    void build_test(uint n_indiv);
	
    void assign_dox_distribution(string dol_distrib,
                                 string doi_distrib,
                                 string doh_distrib);
    
    /**
     * Calculate humoral immunity index for all individuals.
     */
    void assign_immunity_hum();

    /**
     * Calculate cellular immunity index for all individuals.
     */
    void assign_immunity_cell();
    
    void assign_frailty();
    
	// Simulate
	void run();
	
	
	// Set functions
	
    void set_start_time(double t) {_start_time = t;}
    void set_initial_prevalence(uint i0) {_initial_prevalence = i0;}
	void set_current_time(double t) {_current_time = t;}
	void set_world(world w);
	void set_horizon(double h) {_horizon = h;}
	void set_modelParam(modelParam mp) {_modelParam = mp;}
    
    
    /** 
     * Make the inventory of all 'sp_other' social places and record in '_sp_other'.
     */
    void set_sp_other();
    
    void set_sp_other_link(uint k, uint pos_indiv, socialPlace &sp);
	
	// Get functions
	
	world           get_world() {return _world;}
	
	double			get_current_time()   const {return _current_time;}
	vector<double>	get_ts_times()       const {return _ts_times;}
	vector<uint>	get_ts_incidence()   const {return _ts_incidence;}
    
    vector<double>  get_ts_census_sp_time()  const{return _ts_census_sp_time;}
    vector<uint>    get_ts_census_sp_id()    const{return _ts_census_sp_id;}
    vector<string>  get_ts_census_sp_type()  const{return _ts_census_sp_type;}
    vector<uint>    get_ts_census_sp_nS()    const{return _ts_census_sp_nS;}
    vector<uint>    get_ts_census_sp_nE()    const{return _ts_census_sp_nE;}

    uint            get_id_sp_other(uint pos) const {return _sp_other[pos]->get_id_sp();}
    
    vector<uint>    get_track_n_contacts()      const {return _track_n_contacts;}
    vector<uint>    get_track_n_contacts_time() const {return _track_n_contacts_time;}
    vector<uint>    get_track_n_contacts_uid()  const {return _track_n_contacts_uid;}
    
    vector<vector<float> > get_wiw_ages() const {return _wiw_ages;}
    
    vector<vector<double> > get_contactAssort() const {return to_vector_vector(_contactAssort);}
    
	// Time updates
	
	void	time_update(double dt);
	void	update_pop_count();
    void    timeseries_update();
	
	// Migration

	void move_individuals_sched(uint idx_timeslice, double proba);
	void move_individuals(const SPtype sptype, double proba);
	void move_one_individual(uint pos_indiv, ID from, ID to);
    void assign_hospital_to_individuals();
    void change_rnd_sp_other();
    
    socialPlace* pick_rnd_sp_other();
    
    // Book keeping

    void    check_book_keeping();
    
    /**
     * Define all IDs and pointers
     * of tracked individuals
     * in all social places.
     */
    void    define_all_id_tables();
    void    check_schedules_consistency();
	

    // Demographics
    
    void    define_contactAssort();
   
    
	// Epidemic
    
	void	set_disease(const disease& d);
    /**
     * Seed infection in specified socialplaces, with specified initial number of infectious individuals.
     */
	void	seed_infection(vector<ID> id_sp, vector<uint> I0);
	uint	transmission_oneSP(uint k, double dt);
	void	transmission_world(double timeslice);
	uint	prevalence();

    /**
     * Draw the contact rate for an infectious individual
     * in a given social place.
     */
    double          draw_contact_rate(individual* indiv, uint k);
    vector<uint>    draw_n_contacts(uint k,
                                    double dt,
                                    string infectious_type);
    
    vector< vector<uint> > draw_contacted_S(uint k,
                                            vector<uint> n_contacts,
                                            string infectious_type);
    
    vector< vector<uint> > draw_contacted_S_age_constraint(uint k,
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
    void    activate_interventions(ID id_sp, double dt,
                                   float treat_doi_reduc,
                                   float vax_imm_hum_incr,
                                   float vax_imm_cell_incr,
                                   float vax_frail_incr,
                                   float vax_lag);
    
   
    
    
    void    update_immunities();
    bool    at_least_one_vaccination_intervention();
    
    /** Draw the targeted individuals of the ith intervention,
     *  in social place with ID 'id_sp'.
     */
    vector<individual*> draw_targeted_individuals(uint i_intervention,
                                                  ID id_sp,
                                                  double dt);
    
    
    // Census
    
    bool    at_least_one_infected();
    uint    census_total_alive();
    void    update_ts_census_SP();
    
    /**
     * Loop through all SP to get their type and numbre of individuals linked
     */
    dcDataFrame    census_sp();
    
    
	// Exports
	
	dcDataFrame		timeseries();
    
    // Miscelleanous

	uint	population_size();
    void    display_summary_info();
	void	display_split_pop_present();
	void	display_split_pop_linked();
	void	displayInfo_indiv();
	
	void test();
	
	
	
};

double age_contact_elem(double x,double y,double a,double b,double w,double q,double r);


#endif /* defined(__naf__simulation__) */



