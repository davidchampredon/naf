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
	
	// insert row names (as an additional column,
	// don't know how to do differently)
	if(addRowName) {
		rcpplist.push_back(df.get_rowname());
		vector<string> cn = df.get_colname();
		cn.push_back("rowname");
		rcpplist.attr("names") = cn;
	}
	return rcpplist;
}



// [[Rcpp::export]]
List naf_test(List params, List simulParams){
	/// Main test
	
	vector<unsigned int> res;
	
	// Model parameters:
	
	bool debug_mode		= params["debug_mode"];
	
	string dol_distrib	= params["dol_distrib"];
	string doi_distrib	= params["doi_distrib"];
	string doh_distrib	= params["doh_distrib"];
	
	double dol_mean		= params["dol_mean"];
	double doi_mean		= params["doi_mean"];
	double doh_mean		= params["doh_mean"];
	
	double proba_move   = params["proba_move"];
	double contact_rate = params["contact_rate"];
	
	double asymptom_infectiousness_ratio = params["asymptom_infectiousness_ratio"];
	bool homog				= params["homogeneous_contact"];
	double doi_reduc_treat	= params["doi_reduc_treat"];
	double treat_reduc_infect_mean	= params["treat_reduc_infect_mean"];
	

	
	// Simulation parameters:
	
	unsigned int rnd_seed	= simulParams["rnd_seed"];
	unsigned int n_indiv	= simulParams["n_indiv"];
	unsigned int i0			= simulParams["initial_latent"];
	double horizon			= simulParams["horizon"];
	unsigned int nt			= simulParams["nt"];
	unsigned int popexport	= simulParams["popexport"];
	
	string	interv_name		= simulParams["interv_name"];
	string	interv_type		= simulParams["interv_type"];
	string	interv_target	= simulParams["interv_target"];
	double	interv_start	= simulParams["interv_start"];
	double	interv_end  	= simulParams["interv_end"];
	double	interv_cvg_rate  	= simulParams["interv_cvg_rate"];
	double	interv_cvg_max_prop = simulParams["interv_cvg_max_prop"];
	
	
	// Store parameters in a 'modelParam' class:
	
	modelParam MP;
	MP.add_prm_bool("debug_mode", debug_mode);
	MP.add_prm_double("horizon", horizon);
	MP.add_prm_string("dol_distrib", dol_distrib);
	MP.add_prm_string("doi_distrib", doi_distrib);
	MP.add_prm_string("doh_distrib", doh_distrib);
	MP.add_prm_double("dol_mean", dol_mean);
	MP.add_prm_double("doi_mean", doi_mean);
	MP.add_prm_double("doh_mean", doh_mean);
	MP.add_prm_double("proba_move", proba_move);
	MP.add_prm_double("contact_rate", contact_rate);
	MP.add_prm_double("asymptom_infectiousness_ratio", asymptom_infectiousness_ratio);
	MP.add_prm_uint("n_indiv", n_indiv);
	MP.add_prm_bool("homogeneous_contact", homog);
	MP.add_prm_uint("nt", nt);
	MP.add_prm_double ("doi_reduc_treat", doi_reduc_treat);
	MP.add_prm_double ("treat_reduc_infect_mean", treat_reduc_infect_mean);
	
	
	
	cout << "DEBUG: popexport is "<< popexport <<endl;
	
	
	// Define intervention
	intervention interv_test(interv_type,
							 interv_target,
							 interv_name,
							 interv_start,
							 interv_end,
							 interv_cvg_rate,
							 interv_cvg_max_prop);

	
	// Call C++ function
	
	_RANDOM_GENERATOR.seed(rnd_seed);
	
	Simulation sim = test_naf(MP,
							  horizon,
							  n_indiv,
							  i0,
							  interv_test);
	
	// Retrieve all results from simulation:
	// populations:
	dcDataFrame pop_final = sim.get_world()[popexport].export_dcDataFrame();
	// epidemic time series
	dcDataFrame ts = sim.timeseries();
	dcDataFrame ts_census = sim.get_ts_census_by_SP();
	
	// Return R-formatted result:
	return List::create(Named("population_final") = dcDataFrameToRcppList(pop_final,false),
						Named("time_series") = dcDataFrameToRcppList(ts, false),
						Named("time_series_census") = dcDataFrameToRcppList(ts_census, false));
}








// [[Rcpp::export]]
List naf_test_SEIR_vs_ODE(List params, List simulParams){
	/// This runs a SEIR model implemented
	/// in this agent-based model.
	/// The goal is to compare its output to an ODE model.
	
	vector<unsigned int> res;
	
	// Model parameters:
	
	bool debug_mode		= params["debug_mode"];
	double dol_mean		= params["dol_mean"];
	double doi_mean		= params["doi_mean"];
	
	double proba_move   = params["proba_move"];
	double contact_rate = params["contact_rate"];
	bool homog          = params["homogeneous_contact"];
	
	unsigned int rnd_seed	= params["rnd_seed"];
	
	// Simulation parameters:
	unsigned int n_indiv = simulParams["n_indiv"];
	unsigned int i0	     = simulParams["initial_latent"];
	double horizon       = simulParams["horizon"];
	unsigned int nt	     = simulParams["nt"];
	
	// Store parameters in a 'modelParam' class:
	modelParam MP;
	MP.add_prm_bool("debug_mode", debug_mode);
	MP.add_prm_double("dol_mean", dol_mean);
	MP.add_prm_double("doi_mean", doi_mean);
	MP.add_prm_double("proba_move", proba_move);
	MP.add_prm_double("contact_rate", contact_rate);
	MP.add_prm_uint("n_indiv", n_indiv);
	MP.add_prm_bool("homogeneous_contact", homog);
	MP.add_prm_uint("nt", nt);
	
	cout << "DEBUG: the seed is "<<rnd_seed <<endl;
	
	// Set random seed
	_RANDOM_GENERATOR.seed(rnd_seed);
	
	// Call C++ function
	Simulation sim = test_SEIR_vs_ODE(MP,
									  horizon,
									  n_indiv,
									  i0);
	
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
List naf_test_2_sp(List params, List simulParams){
	/// Test mvements between 2 social places
	
	vector<unsigned int> res;
	
	// Model parameters:
	
	bool debug_mode		= params["debug_mode"];
	double dol_mean		= params["dol_mean"];
	double doi_mean		= params["doi_mean"];
	
	double proba_move   = params["proba_move"];
	double contact_rate = params["contact_rate"];
	bool homog          = params["homogeneous_contact"];
	
	unsigned int rnd_seed	= params["rnd_seed"];
	
	// Simulation parameters:
	unsigned int n_indiv = simulParams["n_indiv"];
	unsigned int i0	     = simulParams["initial_latent"];
	double horizon       = simulParams["horizon"];
	unsigned int nt	     = simulParams["nt"];
	
	// Store parameters in a 'modelParam' class:
	modelParam MP;
	MP.add_prm_bool("debug_mode", debug_mode);
	MP.add_prm_double("dol_mean", dol_mean);
	MP.add_prm_double("doi_mean", doi_mean);
	MP.add_prm_double("proba_move", proba_move);
	MP.add_prm_double("contact_rate", contact_rate);
	MP.add_prm_uint("n_indiv", n_indiv);
	MP.add_prm_bool("homogeneous_contact", homog);
	MP.add_prm_uint("nt", nt);
	
	cout << "DEBUG: the seed is "<<rnd_seed <<endl;
	
	// Set random seed
	_RANDOM_GENERATOR.seed(rnd_seed);
	
	// Call C++ function
	Simulation sim = test_move_2_sp(MP,
									horizon,
									n_indiv,
									i0);
	
	// Retrieve all results from simulation:
	// populations:
	dcDataFrame pop_final = sim.get_world()[0].export_dcDataFrame();
	// epidemic time series
	dcDataFrame ts = sim.timeseries();
	dcDataFrame ts_census = sim.get_ts_census_by_SP();
	// Return R-formatted result:
	return List::create(Named("population_final") = dcDataFrameToRcppList(pop_final,false),
						Named("time_series") = dcDataFrameToRcppList(ts, false),
						Named("time_series_census") = dcDataFrameToRcppList(ts_census, false));
}



// [[Rcpp::export]]
List naf_test_hosp(List params, List simulParams){
	/// Test hospitalization
	
	vector<unsigned int> res;
	
	// Model parameters:
	
	bool debug_mode		= params["debug_mode"];
	
	string dol_distrib	= params["dol_distrib"];
	string doi_distrib	= params["doi_distrib"];
	string doh_distrib	= params["doh_distrib"];
	
	double dol_mean		= params["dol_mean"];
	double doi_mean		= params["doi_mean"];
	double doh_mean		= params["doh_mean"];
	
	double proba_move   = params["proba_move"];
	double contact_rate = params["contact_rate"];
	
	double asymptom_infectiousness_ratio = params["asymptom_infectiousness_ratio"];
	bool homog          = params["homogeneous_contact"];
	
	unsigned int rnd_seed	= params["rnd_seed"];
	
	// Simulation parameters:
	unsigned int n_indiv	= simulParams["n_indiv"];
	unsigned int i0		= simulParams["initial_latent"];
	double horizon		= simulParams["horizon"];
	unsigned int nt		= simulParams["nt"];
	
	unsigned int popexport	= simulParams["popexport"];
	
	// Store parameters in a 'modelParam' class:
	modelParam MP;
	MP.add_prm_bool("debug_mode", debug_mode);
	MP.add_prm_double("horizon", horizon);

	MP.add_prm_string("dol_distrib", dol_distrib);
	MP.add_prm_string("doi_distrib", doi_distrib);
	MP.add_prm_string("doh_distrib", doh_distrib);
	
	MP.add_prm_double("dol_mean", dol_mean);
	MP.add_prm_double("doi_mean", doi_mean);
	MP.add_prm_double("doh_mean", doh_mean);
	MP.add_prm_double("proba_move", proba_move);
	MP.add_prm_double("contact_rate", contact_rate);
	MP.add_prm_double("asymptom_infectiousness_ratio", asymptom_infectiousness_ratio);
	MP.add_prm_uint("n_indiv", n_indiv);
	MP.add_prm_bool("homogeneous_contact", homog);
	MP.add_prm_uint("nt", nt);

	
	cout << "DEBUG: popexport is "<< popexport <<endl;
	
	// Set random seed
	_RANDOM_GENERATOR.seed(rnd_seed);
	
	// Call C++ function
	Simulation sim = test_hospitalization(MP,horizon,n_indiv,i0);
	
	// Retrieve all results from simulation:
	// populations:
	dcDataFrame pop_final = sim.get_world()[popexport].export_dcDataFrame();
	// epidemic time series
	dcDataFrame ts = sim.timeseries();
	dcDataFrame ts_census = sim.get_ts_census_by_SP();
	// Return R-formatted result:
	return List::create(Named("population_final") = dcDataFrameToRcppList(pop_final,false),
					  Named("time_series") = dcDataFrameToRcppList(ts, false),
					  Named("time_series_census") = dcDataFrameToRcppList(ts_census, false));
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