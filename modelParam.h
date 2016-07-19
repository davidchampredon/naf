//
//  modelParam.h
//  naf
//
//  Created by David CHAMPREDON on 2016-07-07.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#ifndef __naf__modelParam__
#define __naf__modelParam__

#include <stdio.h>
#include <vector>
#include <string>

#include "dcTools.h"

using namespace std;

class modelParam{
	
	vector<double>		_prm_double_val;
	vector<string>		_prm_double_name;

	vector<unsigned int>	_prm_uint_val;
	vector<string>			_prm_uint_name;
	
	vector<int>			_prm_int_val;
	vector<string>			_prm_int_name;

	vector<bool>			_prm_bool_val;
	vector<string>			_prm_bool_name;

	
	vector< vector<double> >	_prmvec_double_val;
	vector<string>				_prmvec_double_name;

	vector< vector<unsigned int> >	_prmvec_uint_val;
	vector<string>					_prmvec_uint_name;

	vector< vector<int> >		_prmvec_int_val;
	vector<string>				_prmvec_int_name;

public:
	
	modelParam() {}
	
	void add_prm_double(string name, double value)
	{_prm_double_name.push_back(name); _prm_double_val.push_back(value);}
	
	void add_prm_uint(string name, double value)
	{_prm_uint_name.push_back(name); _prm_uint_val.push_back(value);}
	
	void add_prm_int(string name, double value)
	{_prm_int_name.push_back(name); _prm_int_val.push_back(value);}
	
	void add_prm_bool(string name, bool value)
	{_prm_bool_name.push_back(name); _prm_bool_val.push_back(value);}
	
	void add_prmvec_double(string name, vector<double> value)
	{_prmvec_double_name.push_back(name); _prmvec_double_val.push_back(value);}
	
	void add_prmvec_uint(string name, vector<unsigned int> value)
	{_prmvec_uint_name.push_back(name); _prmvec_uint_val.push_back(value);}
	
	void add_prmvec_int(string name, vector< int> value)
	{_prmvec_int_name.push_back(name); _prmvec_int_val.push_back(value);}
	
	
	// Get functions:
	
	double			get_prm_double(string name);
	unsigned int	get_prm_uint(string name);
	int				get_prm_int(string name);
	bool			get_prm_bool(string name);
	
	vector<double>			get_prmvec_double(string name);
	vector<unsigned int>	get_prmvec_uint(string name);
	vector<int>			get_prmvec_int(string name);
};



#endif /* defined(__naf__modelParam__) */
