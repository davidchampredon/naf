//
//  modelParam.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-07.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "modelParam.h"



double modelParam::get_prm_double(string name){
	auto pos = distance(_prm_double_name.begin(), find(_prm_double_name.begin(), _prm_double_name.end(), name));
	stopif(pos >= _prm_double_name.size(), "Parameter '" + name + "' not found!");
	return(_prm_double_val[pos]);
}

unsigned int modelParam::get_prm_uint(string name){
	auto pos = distance(_prm_uint_name.begin(), find(_prm_uint_name.begin(), _prm_uint_name.end(), name));
	stopif(pos >= _prm_uint_name.size(), "Parameter '" + name + "' not found!");
	return(_prm_uint_val[pos]);
}


 int modelParam::get_prm_int(string name){
	auto pos = distance(_prm_int_name.begin(), find(_prm_int_name.begin(), _prm_int_name.end(), name));
	stopif(pos >= _prm_int_name.size(), "Parameter '" + name + "' not found!");
	return(_prm_int_val[pos]);
}


bool modelParam::get_prm_bool(string name){
	auto pos = distance(_prm_bool_name.begin(), find(_prm_bool_name.begin(), _prm_bool_name.end(), name));
	stopif(pos >= _prm_bool_name.size(), "Parameter '" + name + "' not found!");
	return(_prm_bool_val[pos]);
}



vector<double> modelParam::get_prmvec_double(string name){
	auto pos = distance(_prmvec_double_name.begin(), find(_prmvec_double_name.begin(), _prmvec_double_name.end(), name));
	stopif(pos >= _prmvec_double_name.size(), "Parameter '" + name + "' not found!");
	return(_prmvec_double_val[pos]);
}