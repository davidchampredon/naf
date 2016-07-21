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
#include "probaDistribution.h"
#include "dcDataFrame.h"

string SPtype2string(SPtype x);
SPtype int2SPtype(unsigned int i);


class socialPlace: public areaUnit{
    
protected:
    
    SPtype			_type;
    ID				_id_sp;
    unsigned int	_size;
    
    vector<individual> _indiv;
    
    vector<ID>		_linked_indiv_id;  // ID of individuals that are linked to this social place
    
    unsigned int	_prevalence;
    
    unsigned int	_n_S;	// number of susceptible
    unsigned int	_n_E;	// number of individuals in latent stage
    unsigned int	_n_Is;	// number of infectious, symptomatic individuals
    unsigned int	_n_Ia;	// number of infectious, asymptomatic individuals
    unsigned int	_n_R;	// number of individuals who recovered from the disease
    unsigned int	_n_D;	// number of dead individuals
    unsigned int	_n_H;	// number of hospitalized individuals
    
    // DELETE WHEN _indiv_S etc set up (?)
    vector<ID>      _id_S;   // IDs of all susceptible in this social place
    vector<ID>      _id_Is;  // IDs of all sympt. indiv in this social place
    vector<ID>      _id_Ia;  // IDs of all asympt. indiv in this social place
    // -----
    

public:
    
    vector<individual*>     _indiv_S;
    vector<individual*>     _indiv_Is;
    vector<individual*>     _indiv_Ia;
    
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
    
    void set_prevalence(unsigned int p) {_prevalence = p;}
    void set_type(SPtype t){_type = t;}
    void set_n_E(ID n) {_n_E = n;}
    
    void set_id_S(vector<ID> x)  {_id_S  = x;}
    void set_id_Is(vector<ID> x) {_id_Is = x;}
    void set_id_Ia(vector<ID> x) {_id_Ia = x;}
    
    void increase_prevalence() {_prevalence++;}
    
    
    // Get functions:
    
    SPtype			get_type()		const {return _type;}
    ID				get_id_sp()		const {return _id_sp;}
    unsigned long	get_size()		const {return _size;}
    
    unsigned int	get_n_S()		const {return _n_S;}
    unsigned int	get_n_E()		const {return _n_E;}
    unsigned int	get_n_Is()		const {return _n_Is;}
    unsigned int	get_n_Ia()		const {return _n_Ia;}
    unsigned int	get_n_R()		const {return _n_R;}
    
    vector<ID>      get_id_S()      const {return _id_S;}
    vector<ID>      get_id_Is()     const {return _id_Is;}
    vector<ID>      get_id_Ia()     const {return _id_Ia;}
    
    
    unsigned int	get_prevalence()            const {return _prevalence;}
    vector<individual>	get_indiv()             const {return _indiv;}
    individual		get_indiv(unsigned int pos)	const {return _indiv[pos];}
    vector<ID>		get_linked_indiv_id()       const {return _linked_indiv_id;}
    
    individual*      get_mem_indiv(unsigned int i) {return &_indiv[i];}
    
    // Time
    
    void time_update(double dt);
    
    
    // Movements of individuals
    
    ID find_dest(unsigned int pos, unsigned int idx_timeslice);
    ID find_dest_linked(unsigned int pos,
                        unsigned int idx_timeslice,
                        const vector<individual>& indiv_vec);
    
    void add_indiv(individual& indiv);
    void add_indiv(vector<individual>& indiv);
    void remove_indiv(individual& indiv);
    void remove_indiv(unsigned int pos);
    void remove_indiv(vector<unsigned int> posvec);
    void add_linked_indiv(ID id);
    void remove_linked_indiv(ID id);
    ID   n_linked_indiv(){return (ID)(_linked_indiv_id.size());}
    
    
    // Diseases
    
    void		set_disease_to_all_indiv(const disease & d);
    void		acquireDisease(unsigned int pos);
    vector<unsigned int> pick_rnd_susceptibles(unsigned int num);
    vector<ID>	id_infected_bruteforce();
    
    void		update_epidemic_count(const individual& indiv, string move_type);
    
    
    // Book keeping
    
    void        clear_id_S()     {_id_S.clear();}
    void        clear_id_Is()    {_id_Is.clear();}
    void        clear_id_Ia()    {_id_Ia.clear();}
    
    void        add_id_S(ID x)   {_id_S.push_back(x);}
    void        add_id_Is(ID x)  {_id_Is.push_back(x);}
    void        add_id_Ia(ID x)  {_id_Ia.push_back(x);}
    
    void        remove_id_S(ID x) {removeValue(_id_S, x);}
    void        remove_id_Is(ID x) {removeValue(_id_Is, x);}
    void        remove_id_Ia(ID x) {removeValue(_id_Ia, x);}
    
    
    // Census functions:
    // WARNING: brute force counting, hence slow!
    
    unsigned int	census_alive();
    unsigned int	census_infectious();
    unsigned int	census_disease_stage(string stage);
    vector<ID>		census_disease_stage_ID(string stage);
    
    
    // Exports
    dcDataFrame		export_dcDataFrame();
    
    
    // Miscellenaous:
    void	displayInfo();
};


vector<socialPlace> build_world_random(unsigned int n_sp,
                                       vector<areaUnit> auvec);

vector<socialPlace> build_world_simple(vector<SPtype> spt,
                                       vector<unsigned int> n_sp,
                                       vector< probaDistrib<unsigned int> > p_size,
                                       vector<individual>& indiv,
                                       vector<areaUnit> auvec,
                                       unsigned int seed =12345);

void populate_random_with_indiv(vector<socialPlace>& v,
                                unsigned int total_indiv,
                                vector<schedule> sched);

vector<ID> at_least_one_indiv_present(const vector<socialPlace>& x);


unsigned int  choose_SPtype_random(const vector<socialPlace>& sp, SPtype x);


vector<socialPlace> test_world(double sizereduction = 0.001);

void displayPopulationSize(const vector<socialPlace>& sp);





#endif /* defined(__naf__socialPlace__) */
