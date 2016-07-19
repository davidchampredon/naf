//
//  Rwrap.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-15.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include <Rcpp.h>

using namespace Rcpp;

#include "tests.h"



List dcDataFrameToRcppList(dcDataFrame df,
						   bool addRowName){
	
	/// Converts a dcDataFrame to a Rccp list
	/// (helper function)
	
	unsigned long ncol = df.get_colname().size();
	
	// Translate the 'dcDataFrame' into a R list
	// (convert to data frame in R, I don't know how to do it here in Rcpp):
	Rcpp::List rcpplist;
	// each column of the dcDataFrame is a list:
	for(int j=0; j<ncol; j++)
		rcpplist.push_back(df.get_value().extractColumn(j));
	// set the associated column names:
	if(!addRowName) rcpplist.attr("names") = df.get_colname();
	
	// insert row names (as an additional column, don't know how to do differently)
	if(addRowName) {
		rcpplist.push_back(df.get_rowname());
		vector<string> cn = df.get_colname();
		cn.push_back("rowname");
		rcpplist.attr("names") = cn;
	}
	return rcpplist;
}



// [[Rcpp::export]]
List naf_test(List params){
	
	vector<unsigned int> res;
	
	// Unpack parameters:
	
	bool debug_mode		= params["debug_mode"];
	double dol_mean		= params["dol_mean"];
	double doi_mean		= params["doi_mean"];
	double horizon			= params["horizon"];
	unsigned int n_indiv	= params["n_indiv"];
	double proba_move		= params["proba_move"];
	double contact_rate	= params["contact_rate"];
	bool homog	           = params["homogeneous_contact"];
	
	unsigned int rnd_seed	= params["rnd_seed"];
	
	// Store parameters in a 'modelParam' class:
	modelParam MP;
	MP.add_prm_bool("debug_mode", debug_mode);
	MP.add_prm_double("dol_mean", dol_mean);
	MP.add_prm_double("doi_mean", doi_mean);
	MP.add_prm_double("proba_move", proba_move);
	MP.add_prm_double("contact_rate", contact_rate);
	MP.add_prm_uint("n_indiv", n_indiv);
	MP.add_prm_bool("homogeneous_contact", homog);
	
	cout << "DEBUG: the seed is "<<rnd_seed <<endl;
	
	// Set random seed
	_RANDOM_GENERATOR.seed(rnd_seed);
	
	// Call C++ function
	Simulation sim = test_transmission(MP, horizon);
	
	// Retrieve all results from simulation:
	
	// populations:
	dcDataFrame pop_final = sim.get_world()[0].export_dcDataFrame();
	
	// epidemic time series
	dcDataFrame ts = sim.timeseries();
	
	
	// Return R-formatted result:
	return List::create( Named("population_final") = dcDataFrameToRcppList(pop_final,false),
						 Named("time_series") = dcDataFrameToRcppList(ts, false));
}




// [[Rcpp::export]]
List rnd_test(List params){
	
	unsigned long N = 1e1;
	unsigned int rnd_seed = params["rnd_seed"];
	
	_RANDOM_GENERATOR.seed(rnd_seed);
	test_random();
//	std::uniform_real_distribution<double> unif(0.0,1.0);
//	cout << endl << "SEED IS: "<< rnd_seed << endl;
//	for (int i=0; i<N; i++) {
//		double y = unif(_RANDOM_GENERATOR);
//		cout << "TEST-distrib = " << y <<endl;
//	}
	cout << endl;
	
	
	
	_RANDOM_GENERATOR.seed(rnd_seed);
	test_random();
	
	return List::create( Named("res_rnd") = N);
}