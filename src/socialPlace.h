//
//  socialPlace.h
//  naf
//
//  Created by David CHAMPREDON on 2016-07-06.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#ifndef __naf__socialPlace__
#define __naf__socialPlace__

#include <stdio.h>


#include "areaUnit.h"
#include "individual.h"
#include "dcTools.h"
#include "discrete_prob_dist.h"
#include "dcDataFrame.h"

string SPtype2string(SPtype x);
SPtype int2SPtype(uint i);


class socialPlace: public areaUnit{
    
protected:
    
    SPtype  _type;
    ID      _id_sp;
    uint    _size;
    
    vector<individual> _indiv;
    
    vector<ID>		_linked_indiv_id;  // ID of individuals that are linked to this social place
    
    uint	_prevalence;
    
    uint	_n_S;	// number of susceptible
    uint	_n_E;	// number of individuals in latent stage
    uint	_n_Is;	// number of infectious, symptomatic individuals
    uint	_n_Ia;	// number of infectious, asymptomatic individuals
    uint	_n_R;	// number of individuals who recovered from the disease
    uint	_n_D;	// number of dead individuals
    uint	_n_H;	// number of hospitalized individuals
    
    // DELETE WHEN _indiv_S etc set up (?)
    vector<ID>      _id_S;   // IDs of all susceptible in this social place
    vector<ID>      _id_Is;  // IDs of all sympt. indiv in this social place
    vector<ID>      _id_Ia;  // IDs of all asympt. indiv in this social place
    // -----
    

public:
    
    vector<individual*>     _indiv_S;
    vector<individual*>     _indiv_Is;
    vector<individual*>     _indiv_Ia;
    vector<individual*>     _indiv_H;
    vector<individual*>     _indiv_vax;
    
    // Constructors:
    void base_constructor();
    socialPlace();
    socialPlace(ID id, SPtype type);
    socialPlace(ID id_au,
                string name,
                ID id_region,
                string regionName,
                ID id_sp,
                SPtype type);
    socialPlace(areaUnit AU, ID id_sp, SPtype type);
    
    
    // Set functions
    
    void set_prevalence(uint p) {_prevalence = p;}
    void set_type(SPtype t){_type = t;}
    void set_n_E(uint n) {_n_E = n;}
    void set_n_H(uint n) {_n_H = n;}
    void increment_n_H() {_n_H++;}
    void increment_n_D() {_n_D++;}
    void decrement_n_H() {_n_H--;}
    void decrement_n_Is(){_n_Is--;}
    
    void set_id_S(vector<ID> x)  {_id_S  = x;}
    void set_id_Is(vector<ID> x) {_id_Is = x;}
    void set_id_Ia(vector<ID> x) {_id_Ia = x;}
    
    void set_id_sp_hospital(uint pos, socialPlace sp_hospital)
    {_indiv[pos].set_id_sp_hospital(sp_hospital);}
    
    void increase_prevalence() {_prevalence++;}
    
    void set_schedule_indiv(uint pos, schedule sched);
    
    void set_dol_distrib(uint pos, string distrib) {_indiv[pos].set_dol_distrib(distrib);}
    void set_doi_distrib(uint pos, string distrib) {_indiv[pos].set_doi_distrib(distrib);}
    void set_doh_distrib(uint pos, string distrib) {_indiv[pos].set_doh_distrib(distrib);}
    
    /**
     * Sets humoral immunity of the individual in position 'pos' in this social place.
     */
    void set_immunity_hum(uint pos, float val)  {_indiv[pos].set_immunity_hum(val);}
    
    /**
     * Sets cellular immunity of the individual in position 'pos' in this social place.
     */
    void set_immunity_cell(uint pos, float val) {_indiv[pos].set_immunity_cell(val);}
    
    void set_frailty(uint pos, float val)  {_indiv[pos].set_frailty(val);}
    
    void set_sp_other_link(uint pos, socialPlace &sp){_indiv[pos].set_id_sp_other(sp);}
    
    // Get functions:
    
    SPtype			get_type()		const {return _type;}
    ID				get_id_sp()		const {return _id_sp;}
    unsigned long	get_size()		const {return _size;}
    
    uint	get_n_S()		const {return _n_S;}
    uint	get_n_E()		const {return _n_E;}
    uint	get_n_Is()		const {return _n_Is;}
    uint	get_n_Ia()		const {return _n_Ia;}
    uint	get_n_R()		const {return _n_R;}
    uint	get_n_H()		const {return _n_H;}
    uint	get_n_D()		const {return _n_D;}
    
    vector<ID>      get_id_S()      const {return _id_S;}
    vector<ID>      get_id_Is()     const {return _id_Is;}
    vector<ID>      get_id_Ia()     const {return _id_Ia;}
    
    
    uint        get_prevalence()        const {return _prevalence;}
    vector<individual>	get_indiv()     const {return _indiv;}
    individual  get_indiv(uint pos)     const {return _indiv[pos];}
    vector<ID>  get_linked_indiv_id()   const {return _linked_indiv_id;}
    ID          get_linked_indiv_id(uint pos) const {return _linked_indiv_id[pos];}
    
    individual*      get_mem_indiv(uint i) {return &_indiv[i];}
    
    // Time
    
    void time_update(double dt);
    
    
    // Movements of individuals
    
    
    /** Find the ID of the social place the individual is supposed to move to
      * at the timeslice 'idx_timeslice' of the schedule.
      * (individual is in position 'pos' in the vector '_indiv')
     */
    ID find_dest(uint pos, uint idx_timeslice);
    
    ID find_dest_linked(uint pos,
                        uint idx_timeslice,
                        const vector<individual>& indiv_vec);
    
    void add_indiv(individual& indiv);
    
    
    /** Remove an individual given its POSITION
      * in vector '_indiv' in this social place. */
    void remove_indiv(uint pos);
    
    void add_linked_indiv(ID id);
    void remove_linked_indiv(ID id);
    ID   n_linked_indiv() const {return (ID)(_linked_indiv_id.size());}
    
    
    // Diseases
    
    void		set_disease_to_all_indiv(const disease & d);
    void		acquireDisease(uint pos);
    vector<uint> pick_rnd_susceptibles(uint num);
    vector<ID>	id_infected_bruteforce();
    
    void		update_epidemic_count(const individual& indiv, string move_type);
    
    
    // Book keeping
    
    void        clear_id_S()     {_id_S.clear();}
    void        clear_id_Is()    {_id_Is.clear();}
    void        clear_id_Ia()    {_id_Ia.clear();}
    
    void        add_id_S(ID x)   {_id_S.push_back(x);}
    void        add_id_Is(ID x)  {_id_Is.push_back(x);}
    void        add_id_Ia(ID x)  {_id_Ia.push_back(x);}
    
    void        remove_id_S(ID x)  { stopif(_id_S.size()==0, "_id_S is empty and shouldn't be.") ;  removeValue(_id_S, x);}
    void        remove_id_Is(ID x) { stopif(_id_Is.size()==0,"_id_Is is empty and shouldn't be.") ; removeValue(_id_Is, x);}
    void        remove_id_Ia(ID x) { stopif(_id_Ia.size()==0,"_id_Ia is empty and shouldn't be.") ; removeValue(_id_Ia, x);}
    
    uint        find_indiv_pos(ID id);
    uint        find_indiv_X_pos(ID id, string X);
    
    // Census functions:
    // WARNING: brute force counting, hence slow!
    
    uint        census_alive();
    uint        census_alive_age(double age_lo, double age_hi);
    uint        census_infectious();
    uint        census_disease_stage(string stage);
    vector<ID>  census_disease_stage_ID(string stage);
    uint        census_schedule(string name);
    
    
    // Exports
    
    /**
     * Export this social place to a dcDataFrame
     */
    dcDataFrame export_dcDataFrame() const;
    
    
    /**
     * Export this social place to a dcDataFrame. Light version.
     */
    dcDataFrame export_dcDataFrame_light() const;

    // Miscellenaous:
    void	displayInfo();
};



/**
 * Export the vector of social places to a vector of dcDataFrame
 */
vector<dcDataFrame> export_dcDataFrame(const vector<socialPlace> & df);


/**
 * Export the vector of social places to a single dcDataFrame
 */
dcDataFrame export_world(const vector<socialPlace>& x, bool light_output);


vector<socialPlace> build_world_random(uint n_sp,
                                       vector<areaUnit> auvec);

vector<socialPlace> build_world_simple(vector<SPtype> spt,
                                       vector<uint> n_sp,
                                       vector< discrete_prob_dist<uint> > p_size,
                                       vector<individual>& indiv,
                                       vector<areaUnit> auvec,
                                       uint seed =12345);

vector<socialPlace> build_world_simple_2(vector<individual>& indiv,
                                         vector<areaUnit> auvec,
                                         vector<schedule> sched,
                                         uint seed =12345);


void populate_random_with_indiv(vector<socialPlace>& v,
                                uint total_indiv,
                                vector<schedule> sched);

/**
 * Return the IDs of the social places that have
 * at least one indivial present
 */
vector<ID> at_least_one_indiv_present(const vector<socialPlace>& x);


uint  choose_SPtype_random(const vector<socialPlace>& sp, SPtype x);


vector<socialPlace> test_world(double sizereduction = 0.001);

void displayPopulationSize(const vector<socialPlace>& sp);

uint world_size(const vector<socialPlace>& sp);


uint census_schedule(vector<socialPlace> w, string name);

/**
 * Return (in a vector) the ages of all individual in this world.
 */
vector<float> census_ages(vector<socialPlace> w);


#endif /* defined(__naf__socialPlace__) */
