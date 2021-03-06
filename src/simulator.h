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
    
    // Number of individuals treated or vaccinated
    // for each intervention:
    vector<uint>    _n_treated;
    vector<uint>    _n_vaccinated;
    
    
    // Number of individuals by age (i.e. age distribution).
    // example: _n_age[12] = number of indiv whose round(age)=12
    vector<uint> _n_age;
    
    // Calculate the age distribution
    void calc_n_age();
	
    // Age distribution of vaccinated individuals:
    vector<uint>   _n_vaccinated_age;    // number of individuals
    vector<double> _prop_vaccinated_age; // proportion within age
    
    
    
    
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
    vector<uint>         _max_cvg_interv; // Maximum number of indiv receiving associated intervention

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
    
    /** Average frailty at the whole population level. */
    float   _frailty_average;
    
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
    
    /**
     * Calculate frailty index for all individuals.
     */
    void assign_frailty();
    
    /** Calculate average frailty for the whole population. */
    void calc_frailty_average();
    
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
    
    
    /** link social place 'other' whose ID is 'id_sp'
     *  with the individual in _position_ 'pos' currently in kth social place
     */
    void set_sp_other_link(uint k, uint pos_indiv, socialPlace &sp);
	
	// Get functions
	
	world           get_world() {return _world;}
	
	double			get_current_time()   const {return _current_time;}
	vector<double>	get_ts_times()       const {return _ts_times;}
	vector<uint>	get_ts_incidence()   const {return _ts_incidence;}
    vector<uint>	get_ts_n_vaccinated() const {return _ts_n_vaccinated;}
    
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
    
    /** Update epidemiological time series (at each time step)
     */
    void    timeseries_update();
	
	// Migration

    /** Move individuals across social places according to their schedule
    */
	void move_individuals_sched(uint idx_timeslice,
                                double proba,
                                double red_sympt);
	
    
    /** Move the individual in position "pos_indiv" in thevector "_indiv"
      * from one social place to another.
      * (social places are identified by their IDs/position)
     */
	void move_one_individual(uint pos_indiv, ID from, ID to);
    
    void assign_hospital_to_individuals();
    
    /** Change randomly the link with SP_other for all individuals
     */
    void change_rnd_sp_other();
    
    /** Choose randomly a social place of type 'other' */
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

    uint    final_size();
    
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
    
    /** Calculate probability of transmission given contact
      * between an infectious and susceptible individuals.
    */
    double calc_proba_transmission(individual* infectious,
                                            individual* susceptible);
    
    /** Probability to be symptomatic given an individual's immunity and frailty. */
    double calc_proba_symptomatic(float immunity, float frailty,
                                  bool is_vaccinated, double vax_efficacy);
    
    /** Probability to be symptomatic given
      * an individual's frailty.
     */
    double calc_proba_hospitalized(float frailty);

    /** Probability to die at the end of hospitalization period.
     */
    double calc_proba_death(float frailty);
    
    /** Full transmission process:
      * draw number of contacts, identify susceptible contacted,
      * attempts transmission, activate successful attempts.
     */
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
    
    /** Scan all infectious symptomatic individuals
      * and hospitalize them if they were meant to and if it's time to do so.
     */
    void    hospitalize();
    
    /** Scan all infectious hospitalized individuals
     * and discharge them if duration of hospitalization > than the one drawn (_doh_drawn)
     */
    void    discharge_hospital(uint idx_timeslice);
    
    void    death_hospital();

    // Interventions
    
    void    add_intervention(intervention x) {_intervention.push_back(x);}
    
    
    void    count_targeted_by_intervention();
    
    /** Activate all interventions for social place 'id_sp'.
     */
    void    activate_interventions(ID id_sp, double dt,
                                   float treat_doi_reduc,
                                   float vax_imm_hum_incr,
                                   float vax_imm_cell_incr,
                                   float frail_min,
                                   float vax_lag);
    
    
    /** Calculate the earliest intervention start time **/
    double earliest_interv_start();
    
    
   /** Change simulation start time based on interventions timing **/
    void optimize_start_time();
    
    
    /** Update immunity and frailty of vaccinated individuals
     * at each time step. Because vaccine takes some time
     * to reach its full efficacy, both immunity and frailty
     * need to be updated during this period of time.
     */
    void    update_immunities();
    
    
    bool    at_least_one_vaccination_intervention();
    
    
    /** Helper function for the "priority_age_fraity" intervention.
     */
    vector<individual*> helper_priority_vax(uint i_intervention,
                                            ID id_sp,
                                            double dt,
                                            float age_young,
                                            float age_old);
    
    /** Draw the targeted individuals of the ith intervention,
     *  in social place with ID 'id_sp'.
     */
    vector<individual*> draw_targeted_individuals(uint i_intervention,
                                                  ID id_sp,
                                                  double dt);
    
    /** retrieve the _first_ vaccine efficacy encounter in all interventions. */
    double retrieve_vaccine_efficacy();
    
    /** Calculate the proportion of vaccinated indiv
     *  of a specified age.
     */
    float calc_cumvax_prop(string rngAge, int age);
    
    
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
    
    /** Size distributions of all SP for each type.
     *  result[i]: vector of all size values for SP type 'i'.
     */
    vector< vector<unsigned long> > sp_size_distribution();

    void    display_summary_info();
	void	display_split_pop_present();
	void	display_split_pop_linked();
	void	displayInfo_indiv();
	
	void test();
	
	
	
};

/** retrieve the minimum _MEAN_ frailty. */
double minimum_frailty(world w, float f0, float agepivot, float slope1, float slope2);
double age_contact_elem(double x,double y,double a,double b,double w,double q,double r);



#endif /* defined(__naf__simulation__) */



