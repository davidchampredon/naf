//
//  disease.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-15.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "disease.h"



void disease::base_constructor(){
	/// Base constructor in all genuine constructors
	_name = "_UNDEFINED_DISEASE_NAME";
	_dol_mean = 0.0;
	_doi_mean =	0.0;
    _doh_mean = 0.0;
}


disease::disease(string name){
	base_constructor();
	_name = name;
}

disease::disease(string name,
                 float dol_mean,
                 float doi_mean){
    base_constructor();
    _name = name;
    _dol_mean = dol_mean;
    _doi_mean = doi_mean;
}


disease::disease(string name,
                 float dol_mean,
                 float doi_mean,
                 float doh_mean){
	base_constructor();
	_name = name;
	_dol_mean = dol_mean;
	_doi_mean = doi_mean;
    _doh_mean = doi_mean;
}

