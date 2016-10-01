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


List to_list_vector(const vector<dcDataFrame> & df,
                    bool addRowName){
    /// Converts a vector of dcDataFrame to a Rccp list
    /// (helper function)
    
    unsigned long n = df.size();
    
    // Translate the 'dcDataFrame' into a R list
    // (convert to data frame in R, I don't know how to do it here in Rcpp):
    Rcpp::List rcpplist;
    
    for(unsigned long i=0; i<n; i++)
        rcpplist.push_back(to_list(df[i],addRowName));
    
    return rcpplist;
}




void set_parameter(modelParam &MP, string prm_name, string prm_type, List params){
    
    bool found = false;
    
    if (prm_type=="double") {
        double val = params[prm_name];
        MP.add_prm_double(prm_name, val);
        found = true;
    }
    
    else if (prm_type=="bool") {
        bool val = params[prm_name];
        MP.add_prm_bool(prm_name, val);
        found = true;
    }
    
    else if (prm_type=="string") {
        string val = params[prm_name];
        MP.add_prm_string(prm_name, val);
        found = true;
    }
    
    else if (prm_type=="uint") {
        uint val = params[prm_name];
        MP.add_prm_uint(prm_name, val);
        found = true;
    }
    if(!found)
        cerr << "Parameter type " << prm_type << " not found for parameter "<<prm_name <<endl;
}



// ====================================================================================
// ====================================================================================
// = = = =  M A I N    E X P O R T  = = = =
// ====================================================================================
// ====================================================================================


// [[Rcpp::export]]
List naf_run(List params,
             List simulParams,
             List intervParams,
             List worldParams,
             List scheduleParams){
    
    try{
        /// Main test
        
        cout << " Entering Rwrap... " << endl;
        
        vector<unsigned int> res;
        
        // ==== Model parameters ====
        
        modelParam MP;
        
        set_parameter(MP, "debug_mode", "bool", params);
        
        set_parameter(MP, "homogeneous_contact", "bool", params);
        
        set_parameter(MP, "dol_distrib", "string", params);
        set_parameter(MP, "doi_distrib", "string", params);
        set_parameter(MP, "doh_distrib", "string", params);
        
        set_parameter(MP, "dol_mean", "double", params);
        set_parameter(MP, "doi_mean", "double", params);
        set_parameter(MP, "doh_mean", "double", params);
        
        set_parameter(MP, "dol_var", "double", params);
        set_parameter(MP, "doi_var", "double", params);
        set_parameter(MP, "doh_var", "double", params);
        
        set_parameter(MP, "proba_move", "double", params);
        set_parameter(MP, "proba_change_sp_other", "double", params);
        
        set_parameter(MP, "contact_rate_mean",             "double", params);
        set_parameter(MP, "contact_rate_stddev",           "double", params);
        set_parameter(MP, "contact_ratio_age_1_10",        "double", params);
        set_parameter(MP, "contact_ratio_age_10_16",       "double", params);
        set_parameter(MP, "contact_ratio_age_over_65",     "double", params);
        set_parameter(MP, "contact_ratio_sp_household",    "double", params);
        set_parameter(MP, "contact_ratio_sp_pubTransport", "double", params);
		
		set_parameter(MP, "contactAssort_lambda", "double", params);
        
        set_parameter(MP, "asymptom_infectiousness_ratio", "double", params);
        
        set_parameter(MP, "treat_doi_reduc", "double", params);
        set_parameter(MP, "treat_reduc_infect_mean", "double", params);
        
        set_parameter(MP, "vax_imm_incr",          "double", params);
        set_parameter(MP, "vax_frail_incr",        "double", params);
        set_parameter(MP, "vax_lag_full_efficacy", "double", params);
        
        set_parameter(MP, "proba_death_prm_1", "double", params);
        set_parameter(MP, "proba_death_prm_2", "double", params);
        set_parameter(MP, "proba_death_prm_3", "double", params);
        
        set_parameter(MP, "frailty_0",         "double", params);
        set_parameter(MP, "frailty_min",       "double", params);
        set_parameter(MP, "frailty_agemin",    "double", params);
        set_parameter(MP, "frailty_agepivot",  "double", params);
        set_parameter(MP, "frailty_pivot",     "double", params);
        set_parameter(MP, "frailty_powerChild","double", params);
        set_parameter(MP, "frailty_sd",        "double", params);
		
		set_parameter(MP, "imm_hum_baseline",  "double", params);
		set_parameter(MP, "imm_hum_agezero",   "double", params);
		set_parameter(MP, "imm_hum_p",         "double", params);
        
        
        // ==== Simulator parameters ====
        
        unsigned int rnd_seed   = simulParams["rnd_seed"];
        unsigned int i0         = simulParams["initial_latent"];
        double horizon          = simulParams["horizon"];
        unsigned int popexport  = simulParams["popexport"];
        double start_time       = simulParams["start_time"];
        
        
        // === Interventions ===
        
        vector<intervention> interv_vec;
        
        unsigned long n_interv = intervParams.size();
        
        cout << "Loading "<< n_interv <<" intervention scenarios ... ";
        
        
        for (unsigned long i=0; i<n_interv; i++) {
            
            Rcpp::List tmp_interv = intervParams[i];
            
            string	interv_name			= tmp_interv["interv_name"];
            string	interv_type			= tmp_interv["interv_type"];
            string	interv_target		= tmp_interv["interv_target"];
            double	interv_start		= tmp_interv["interv_start"];
            double	interv_end			= tmp_interv["interv_end"];
            double	interv_cvg_rate		= tmp_interv["interv_cvg_rate"];
            double	interv_cvg_max_prop = tmp_interv["interv_cvg_max_prop"];

            intervention interv(interv_type,
                                interv_target,
                                interv_name,
                                interv_start,
                                interv_end,
                                interv_cvg_rate,
                                interv_cvg_max_prop);

            // Add to the vector holding all interventions:
            interv_vec.push_back(interv);
            
            if(i>0) cout << ","; cout<< i+1;
        }
        
        cout << endl << "Interventions loaded. " << endl;
        
        
        // === World population parameters ===
        
        vector<ID> id_au        = worldParams["id_au"];
        vector<string> name_au  = worldParams["name_au"];
        ID id_region            = worldParams["id_region"];
        string regionName       = worldParams["regionName"];
        
        
        vector<areaUnit> auvec = create_area_unit(id_au,
                                                  name_au,
                                                  id_region,
                                                  regionName);
        
        vector<uint> n_hh       = worldParams["n_hh"];
        vector<uint> n_wrk      = worldParams["n_wrk"];
        vector<uint> n_pubt     = worldParams["n_pubt"];
        vector<uint> n_school   = worldParams["n_school"];
        vector<uint> n_hosp     = worldParams["n_hosp"];
        vector<uint> n_other    = worldParams["n_other"];
        
        
        // Social place sizes distributions:
        
        vector<discrete_prob_dist<uint> > D_size_hh_vec;
        vector<discrete_prob_dist<uint> > D_size_wrk_vec;
        vector<discrete_prob_dist<uint> > D_size_pubt_vec;
        vector<discrete_prob_dist<uint> > D_size_school_vec;
        vector<discrete_prob_dist<uint> > D_size_hosp_vec;
        vector<discrete_prob_dist<uint> > D_size_other_vec;
        
        unsigned long n_au = id_au.size();
        
        List tmp_szdst;
        
        for(unsigned int i=0; i<n_au; i++){
            

            tmp_szdst = worldParams["hh_size"];
            vector<uint> hh_size            = tmp_szdst[i];

            tmp_szdst = worldParams["hh_size_proba"];
            vector<double> hh_size_proba    = tmp_szdst[i];
            
            tmp_szdst = worldParams["wrk_size"];
            vector<uint> wrk_size           = tmp_szdst[i];
            
            tmp_szdst = worldParams["wrk_size_proba"];
            vector<double> wrk_size_proba   = tmp_szdst[i];
            
            tmp_szdst = worldParams["pubt_size"];
            vector<uint> pubt_size          = tmp_szdst[i];
            
            tmp_szdst = worldParams["pubt_size_proba"];
            vector<double> pubt_size_proba  = tmp_szdst[i];
            
            tmp_szdst = worldParams["school_size"];
            vector<uint> school_size        = tmp_szdst[i];
            
            tmp_szdst = worldParams["school_size_proba"];
            vector<double> school_size_proba= tmp_szdst[i];
            
            tmp_szdst = worldParams["hosp_size"];
            vector<uint> hosp_size          = tmp_szdst[i];
            
            tmp_szdst = worldParams["hosp_size_proba"];
            vector<double> hosp_size_proba  = tmp_szdst[i];
            
            tmp_szdst = worldParams["other_size"];
            vector<uint> other_size         = tmp_szdst[i];
            
            tmp_szdst = worldParams["other_size_proba"];
            vector<double> other_size_proba = tmp_szdst[i];

            discrete_prob_dist<uint> D_size_hh(hh_size, hh_size_proba);
            discrete_prob_dist<uint> D_size_wrk(wrk_size, wrk_size_proba);
            discrete_prob_dist<uint> D_size_pubt(pubt_size, pubt_size_proba);
            discrete_prob_dist<uint> D_size_school(school_size, school_size_proba);
            discrete_prob_dist<uint> D_size_hosp(hosp_size, hosp_size_proba);
            discrete_prob_dist<uint> D_size_other(other_size, other_size_proba);
            
            D_size_hh_vec.push_back(D_size_hh);
            D_size_wrk_vec.push_back(D_size_wrk);
            D_size_pubt_vec.push_back(D_size_pubt);
            D_size_school_vec.push_back(D_size_school);
            D_size_hosp_vec.push_back(D_size_hosp);
            D_size_other_vec.push_back(D_size_other);
        }
        
        // Within households age conditional distributions:
		
		// Retrieve the maximum houshold size across all AU:
		vector<uint> hms = worldParams["max_hh_size"];
		uint hh_mxsz = *max_element(hms.begin(), hms.end());
		
		vector< vector<discrete_prob_dist<uint> > > pr_age_hh;
		
		
		for(uint h=0; h<hh_mxsz; h++){
			
			vector<discrete_prob_dist<uint> > tmp;
			
			for(uint i=0; i<=h; i++){
				string n_v = "pr_age_hh_"+ to_string(h)+to_string(i)+"_val";
				string n_p = "pr_age_hh_"+ to_string(h)+to_string(i)+"_proba";
				
				vector<uint>   pr_age_hh_val   = worldParams[n_v];
				vector<double> pr_age_hh_proba = worldParams[n_p];
				
				discrete_prob_dist<uint> dist(pr_age_hh_val,
											  pr_age_hh_proba);
				
				tmp.push_back(dist);
			}
			pr_age_hh.push_back(tmp);
		}
		
		cout << "all parameters loaded." << endl;
		
        // === Schedules definition ===
        
        // Define the time slices
        // must be same for all schedules and must sum up to 1.0
        vector<double> timeslice           = scheduleParams["timeslice"];
		vector<string> sched_names         = scheduleParams["sched_names"];
		vector<vector<string> > sched_desc = scheduleParams["sched_desc"];
		
		vector<schedule> sched = build_all_schedules(sched_desc, sched_names, timeslice);
		
		
		// =========================
        // === Call C++ function ===
		// =========================
		
        _RANDOM_GENERATOR.seed(rnd_seed);
        
        Simulator sim = run_test(auvec,
                                 D_size_hh_vec,
                                 D_size_wrk_vec,
                                 D_size_pubt_vec,
                                 D_size_school_vec,
                                 D_size_hosp_vec,
                                 D_size_other_vec,
                                 pr_age_hh,
                                 n_hh ,
                                 n_wrk,
                                 n_pubt ,
                                 n_school,
                                 n_hosp,
                                 n_other,
                                 sched ,
                                 MP,
                                 start_time,
                                 horizon,
                                 i0,
                                 interv_vec);
        
        // Retrieve all results from simulation:
        // populations:
        dcDataFrame pop_final = sim.get_world()[popexport].export_dcDataFrame();
        // epidemic time series
        dcDataFrame ts = sim.timeseries();
        
        // Time series of census for every social places:
        Rcpp::List census_sp;
        census_sp.push_back(sim.get_ts_census_sp_time());
        census_sp.push_back(sim.get_ts_census_sp_id());
        census_sp.push_back(sim.get_ts_census_sp_type());
        census_sp.push_back(sim.get_ts_census_sp_nS());
        census_sp.push_back(sim.get_ts_census_sp_nE());
        vector<string> tmp_names = {"time","id_sp","type","nS","nE"};
        census_sp.attr("names") = tmp_names;
        
        // Number of contacts tracker
        Rcpp::List track_n_contacts;
        track_n_contacts.push_back(sim.get_track_n_contacts_time());
        track_n_contacts.push_back(sim.get_track_n_contacts_uid());
        track_n_contacts.push_back(sim.get_track_n_contacts());
        tmp_names = {"time","uid","nContacts"};
        track_n_contacts.attr("names") = tmp_names;
        
        
        vector<dcDataFrame> W = export_dcDataFrame(sim.get_world());
        
        // Return R-formatted result:
        return List::create(Named("population_final") = to_list(pop_final,false),
                            Named("world")            = to_list_vector(W, false),
                            Named("time_series")      = to_list(ts, false),
                            Named("time_series_sp")   = census_sp,
                            Named("track_n_contacts") = track_n_contacts,
                            Named("wiw_ages")         = sim.get_wiw_ages(),
							Named("contactAssort")    = sim.get_contactAssort()
                            );
    }
    catch (...){
        ::Rf_error(">>>> C++ exception (unknown reason) <<<<");
        return NULL;
    }
}





// ====================================================================================
// ====================================================================================
// = = = =  T E S T S  = = = =
// ====================================================================================
// ====================================================================================




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
    
    
    // Simulator parameters:
    
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
    
    Simulator sim = test_naf(MP,
                             horizon,
                             n_indiv,
                             i0,
                             interv_test);
    
    // Retrieve all results from simulation:
    // populations:
    dcDataFrame pop_final = sim.get_world()[popexport].export_dcDataFrame();
    // epidemic time series
    dcDataFrame ts = sim.timeseries();
    
    // Return R-formatted result:
    return List::create(Named("population_final") = to_list(pop_final,false),
                        Named("time_series") = to_list(ts, false),
                        Named("time_series_census") = 0);
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
    
    // Simulator parameters:
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
    Simulator sim = test_SEIR_vs_ODE(MP,
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
    
    // Simulator parameters:
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
    Simulator sim = test_move_2_sp(MP,
                                   horizon,
                                   n_indiv,
                                   i0);
    
    // Retrieve all results from simulation:
    // populations:
    dcDataFrame pop_final = sim.get_world()[0].export_dcDataFrame();
    // epidemic time series
    dcDataFrame ts = sim.timeseries();
    
    // Return R-formatted result:
    return List::create(Named("population_final") = to_list(pop_final,false),
                        Named("time_series") = to_list(ts, false),
                        Named("time_series_census") = 0);
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
    
    // Simulator parameters:
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
    Simulator sim = test_hospitalization(MP,horizon,n_indiv,i0);
    
    // Retrieve all results from simulation:
    // populations:
    dcDataFrame pop_final = sim.get_world()[popexport].export_dcDataFrame();
    // epidemic time series
    dcDataFrame ts = sim.timeseries();
    
    // Return R-formatted result:
    return List::create(Named("population_final") = to_list(pop_final,false),
                        Named("time_series") = to_list(ts, false),
                        Named("time_series_census") = 0);
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