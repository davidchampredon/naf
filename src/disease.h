//
//  disease.h
//  naf
//
//  Created by David CHAMPREDON on 2016-07-15.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#ifndef __naf__disease__
#define __naf__disease__

#include <stdio.h>
#include <string>

using namespace std;

class disease{
	
	string	_name;
	
	float	_dol_mean;		// mean duration of latency
	float	_doi_mean;		// mean duration of infectiousness
	float	_doh_mean;		// mean duration of hospitalization
    
    float   _dol_var;       // variance for duration of latency
    float   _doi_var;       // variance for duration of infectiousness
    float   _doh_var;       // variance for duration of hospitalization
    
public:
	
	disease(){base_constructor();}
	
	void	base_constructor();
	
	disease(string name);
    disease(string name,
            float dol_mean,
            float doi_mean);
	disease(string name,
            float dol_mean,
            float doi_mean,
            float doh_mean);
    
    disease(string name,
            float dol_mean,
            float doi_mean,
            float doh_mean,
            float dol_var,
            float doi_var,
            float doh_var);
	
	// Get functions:
	string	get_name()      const {return _name;}
	float	get_dol_mean()  const {return _dol_mean;}
	float	get_doi_mean()  const {return _doi_mean;}
    float	get_doh_mean()  const {return _doh_mean;}
    float	get_dol_var()   const {return _dol_var;}
    float	get_doi_var()   const {return _doi_var;}
    float	get_doh_var()   const {return _doh_var;}
    
};


#endif /* defined(__naf__disease__) */
