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



List to_list(const dcDataFrame & df,
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
    
    double doi_reduc_treat			= params["treat_doi_reduc"];
    double treat_reduc_infect_mean	= params["treat_reduc_infect_mean"];
    
    double vax_imm_incr         = params["vax_imm_incr"];
    double vax_frail_incr       = params["vax_frail_incr"];
    double vax_lag_full_efficacy= params["vax_lag_full_efficacy"];
    
    double proba_death_prm_1    = params["proba_death_prm_1"];
    double proba_death_prm_2    = params["proba_death_prm_2"];
    double proba_death_prm_3    = params["proba_death_prm_3"];
    
    
    // Simulation parameters:
    
    unsigned int rnd_seed   = simulParams["rnd_seed"];
    unsigned int n_indiv    = simulParams["n_indiv"];
    unsigned int i0         = simulParams["initial_latent"];
    double horizon          = simulParams["horizon"];
    unsigned int nt         = simulParams["nt"];
    unsigned int popexport  = simulParams["popexport"];
    
    string	interv_name     = simulParams["interv_name"];
    string	interv_type     = simulParams["interv_type"];
    string	interv_target   = simulParams["interv_target"];
    double	interv_start    = simulParams["interv_start"];
    double	interv_end      = simulParams["interv_end"];
    double	interv_cvg_rate = simulParams["interv_cvg_rate"];
    double	interv_cvg_max_prop = simulParams["interv_cvg_max_prop"];
    
    
    // Store parameters in a 'modelParam' class:
    
    modelParam MP;
    
    MP.add_prm_bool     ("debug_mode", debug_mode);
    MP.add_prm_double   ("horizon", horizon);
    MP.add_prm_string   ("dol_distrib", dol_distrib);
    MP.add_prm_string   ("doi_distrib", doi_distrib);
    MP.add_prm_string   ("doh_distrib", doh_distrib);
    MP.add_prm_double   ("dol_mean", dol_mean);
    MP.add_prm_double   ("doi_mean", doi_mean);
    MP.add_prm_double   ("doh_mean", doh_mean);
    MP.add_prm_double   ("proba_move", proba_move);
    
    MP.add_prm_double   ("proba_death_prm_1", proba_death_prm_1);
    MP.add_prm_double   ("proba_death_prm_2", proba_death_prm_2);
    MP.add_prm_double   ("proba_death_prm_3", proba_death_prm_3);
    
    MP.add_prm_double   ("contact_rate", contact_rate);
    MP.add_prm_double   ("asymptom_infectiousness_ratio", asymptom_infectiousness_ratio);
    MP.add_prm_uint     ("n_indiv", n_indiv);
    MP.add_prm_bool     ("homogeneous_contact", homog);
    MP.add_prm_uint     ("nt", nt);
    MP.add_prm_double   ("treat_doi_reduc", doi_reduc_treat);
    MP.add_prm_double   ("treat_reduc_infect_mean", treat_reduc_infect_mean);
    
    MP.add_prm_double   ("vax_imm_incr", vax_imm_incr);
    MP.add_prm_double   ("vax_frail_incr", vax_frail_incr);
    MP.add_prm_double   ("vax_lag_full_efficacy", vax_lag_full_efficacy);
    
    
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
    return List::create(Named("population_final") = to_list(pop_final,false),
                        Named("time_series") = to_list(ts, false),
                        Named("time_series_census") = to_list(ts_census, false));
}


// [[Rcpp::export]]
List build_world(List params){
	
	
	vector<ID> id_au = params["id_au"];
	vector<string> name_au = params["name_au"];

	ID id_region = 0;
	string regionName = "RegionOne";
	vector<areaUnit> auvec = create_area_unit(id_au, name_au, id_region, regionName);
	
	// Households sizes
	vector<uint> hh_size			= params["hh_size"];
	vector<double> hh_size_proba	= params["hh_size_proba"];
	discrete_prob_dist<uint> D_size_hh(hh_size, hh_size_proba);
	D_size_hh.display();
	
	// Age distribution inside households
	vector< vector<discrete_prob_dist<uint> > > pr_age_hh;
	
	discrete_prob_dist<uint> pr_age_hh_00(params["age_hh_00"],params["age_hh_00_proba"]);
	discrete_prob_dist<uint> pr_age_hh_10(params["age_hh_10"],params["age_hh_10_proba"]);
	discrete_prob_dist<uint> pr_age_hh_11(params["age_hh_11"],params["age_hh_11_proba"]);
	discrete_prob_dist<uint> pr_age_hh_20(params["age_hh_20"],params["age_hh_20_proba"]);
	discrete_prob_dist<uint> pr_age_hh_21(params["age_hh_21"],params["age_hh_21_proba"]);
	discrete_prob_dist<uint> pr_age_hh_22(params["age_hh_22"],params["age_hh_22_proba"]);
	
	vector<discrete_prob_dist<uint> > tmp;
	
	tmp.push_back(pr_age_hh_00);
	pr_age_hh.push_back(tmp);
	
	tmp.clear();
	tmp = {pr_age_hh_10,pr_age_hh_11};
	pr_age_hh.push_back(tmp);
	
	tmp.clear();
	tmp = {pr_age_hh_20,pr_age_hh_21,pr_age_hh_22};
	pr_age_hh.push_back(tmp);
	
	
	// Workplace sizes
	vector<uint> wrk_size			= params["wrk_size"];
	vector<double> wrk_size_proba	= params["wrk_size_proba"];
	discrete_prob_dist<uint> D_size_wrk(wrk_size, wrk_size_proba);
	D_size_wrk.display();
	
	// Public transport sizes
	vector<uint> pubt_size			= params["pubt_size"];
	vector<double> pubt_size_proba	= params["pubt_size_proba"];
	discrete_prob_dist<uint> D_size_pubt(pubt_size, pubt_size_proba);
	D_size_pubt.display();
	
	// School sizes
	vector<uint> school_size			= params["school_size"];
	vector<double> school_size_proba	= params["school_size_proba"];
	discrete_prob_dist<uint> D_size_school(school_size, school_size_proba);
	D_size_school.display();
	
 
	
	// === Create world ===
	
	// number of each social place type
	vector<uint> n_hh       = params["n_hh"];
	vector<uint> n_wrk      = params["n_wrk"];
	vector<uint> n_pubt     = params["n_pubt"];
	vector<uint> n_school   = params["n_school"];
	
	vector<socialPlace> W = build_world(auvec,
									D_size_hh,      // Households sizes
									pr_age_hh,      // Age distribution inside households
									D_size_wrk,
									D_size_pubt,
									D_size_school,
									n_hh,           // number of households
									n_wrk,
									n_pubt,
									n_school);
	
	//DEBUG
	cout << "DEBUG: W.size = "<< W.size() << endl;
	
	vector<dcDataFrame> df = export_dcDataFrame(W);
	
	Rcpp::List rcpplist;

	for(unsigned int i=0; i< df.size(); i++)
		rcpplist.push_back( to_list(df[i], false) );
	
	return rcpplist;
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
    unsigned int n_indiv    = simulParams["n_indiv"];
    unsigned int i0         = simulParams["initial_latent"];
    double horizon          = simulParams["horizon"];
    unsigned int nt         = simulParams["nt"];
    
    // Store parameters in a 'modelParam' class:
    modelParam MP;
    MP.add_prm_bool     ("debug_mode", debug_mode);
    MP.add_prm_double   ("dol_mean", dol_mean);
    MP.add_prm_double   ("doi_mean", doi_mean);
    MP.add_prm_double   ("proba_move", proba_move);
    MP.add_prm_double   ("contact_rate", contact_rate);
    MP.add_prm_uint     ("n_indiv", n_indiv);
    MP.add_prm_bool     ("homogeneous_contact", homog);
    MP.add_prm_uint     ("nt", nt);
    
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
    return List::create( Named("population_final") = to_list(pop_final,false),
                        Named("time_series") = to_list(ts, false));
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
    return List::create(Named("population_final") = to_list(pop_final,false),
                        Named("time_series") = to_list(ts, false),
                        Named("time_series_census") = to_list(ts_census, false));
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
    return List::create(Named("population_final") = to_list(pop_final,false),
                        Named("time_series") = to_list(ts, false),
                        Named("time_series_census") = to_list(ts_census, false));
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