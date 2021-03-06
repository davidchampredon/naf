//
//  Rwrap.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-15.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include <Rcpp.h>

using namespace Rcpp;

#include "run.h"



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


void set_parameter(modelParam &MP, string prm_name, string prm_type, List params)
{
	bool debug = false;
	if (debug) cout << "setting parameter "<< prm_name << "... ";
	
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
	
	if (debug) cout << "done."<<endl;
}

/** 
 * Read a discrete distribution.
 * (Predefined format mandatory!)
 */
void read_discr_distr(List prm,
					  string name,
					  vector<uint> &val,
					  vector<double> &proba,
					  uint n_au,
					  uint i){

	string namep = name + "_proba";
	
	if(n_au==1){
		vector<uint>   v = prm[name];
		vector<double> p = prm[namep];
		val   = v;
		proba = p;
	}
	else{
		List tmp   = prm[name];
		List tmp_p = prm[namep];
		vector<uint>   v = tmp[i];
		vector<double> p = tmp_p[i];
		val   = v;
		proba = p;
	}
}



// ====================================================================================
// ====================================================================================
// = = = =  M A I N    E X P O R T  = = = =
// ====================================================================================
// ====================================================================================



// [[Rcpp::export]]
List naf_run_det(List params,
				 List simulParams,
				 List intervParams,
				 List worldParams,
				 List scheduleParams,
				 unsigned int rnd_seed){
	
	try{
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
		
		cout << "--> Durations loaded." << endl;
		
		set_parameter(MP, "proba_hosp",             "double", params);
		
		set_parameter(MP, "proba_move",             "double", params);
		set_parameter(MP, "proba_move_reduc_sympt", "double", params);
		set_parameter(MP, "proba_change_sp_other",  "double", params);
		
		set_parameter(MP, "proba_death_min",      "double", params);
		set_parameter(MP, "proba_death_max",      "double", params);
		set_parameter(MP, "proba_death_frailCvx", "double", params);
		set_parameter(MP, "proba_death_slope",    "double", params);
		set_parameter(MP, "proba_death_mult",    "double", params);

		cout << "--> Probabilities loaded." << endl;
		
		set_parameter(MP, "contact_rate_distrib",          "string", params);
		set_parameter(MP, "contact_rate_mean",             "double", params);
		set_parameter(MP, "contact_rate_stddev",           "double", params);
		set_parameter(MP, "contact_rate_CV",               "double", params);
		
		set_parameter(MP, "contact_ratio_age_1_10",        "double", params);
		set_parameter(MP, "contact_ratio_age_10_16",       "double", params);
		set_parameter(MP, "contact_ratio_age_16_25",       "double", params);
		set_parameter(MP, "contact_ratio_age_25_40",       "double", params);
		set_parameter(MP, "contact_ratio_age_over_65",     "double", params);
		set_parameter(MP, "contact_ratio_sp_household",    "double", params);
		set_parameter(MP, "contact_ratio_sp_pubTransport", "double", params);
		set_parameter(MP, "contact_ratio_sp_school",       "double", params);
		
		set_parameter(MP, "contactAssort_lambda", "double", params);
		
		cout << "--> Contact params loaded." << endl;
		
		set_parameter(MP, "mult_proba_symptomatic",        "double", params);
		set_parameter(MP, "asymptom_infectiousness_ratio", "double", params);
		
		cout << "--> Asymptomatic params loaded." << endl;
		
		set_parameter(MP, "treat_doi_reduc", "double", params);
		set_parameter(MP, "treat_reduc_infect_mean", "double", params);
		
		set_parameter(MP, "vax_imm_hum_incr",	   "double", params);
		set_parameter(MP, "vax_imm_cell_incr",     "double", params);
		set_parameter(MP, "vax_lag_full_efficacy", "double", params);
		
		cout << "--> Treatment & vax params loaded." << endl;
		
		set_parameter(MP, "frailty_0",         "double", params);
		set_parameter(MP, "frailty_agepivot",  "double", params);
		set_parameter(MP, "frailty_slope1",    "double", params);
		set_parameter(MP, "frailty_slope2",    "double", params);
		set_parameter(MP, "frailty_sd",        "double", params);
		
		set_parameter(MP, "imm_hum_baseline",  "double", params);
		set_parameter(MP, "imm_hum_agezero",   "double", params);
		set_parameter(MP, "imm_hum_p",         "double", params);
		
		set_parameter(MP, "imm_cell_max",	"double", params);
		set_parameter(MP, "imm_cell_slope",	"double", params);
		set_parameter(MP, "imm_cell_pivot", "double", params);
		
		cout << "--> Frailty & immunities loaded." << endl;
		
		
		// ==== Simulator parameters ====
		
		bool build_world_only   = simulParams["build_world_only"];
		unsigned int i0         = simulParams["initial_latent"];
		double horizon          = simulParams["horizon"];
		unsigned int popexport  = simulParams["popexport"];
		double start_time       = simulParams["start_time"];
		bool light_output       = simulParams["light_output"];
		bool superlight_output  = simulParams["superlight_output"];
		
		cout << "--> Simulation params loaded." << endl;
		
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
			double	interv_efficacy     = tmp_interv["interv_efficacy"];
			double	interv_cvg_age_mult = tmp_interv["interv_cvg_age_mult"];
			
			intervention interv(interv_type,
								interv_target,
								interv_name,
								interv_start,
								interv_end,
								interv_cvg_rate,
								interv_cvg_max_prop,
								interv_efficacy,
								interv_cvg_age_mult);
			
			// Add to the vector holding all interventions:
			interv_vec.push_back(interv);
			
			if(i>0) cout << ","; cout<< i+1;
		}
		
		cout << endl << "Interventions loaded. " << endl;
		
		
		// === World population parameters ===
		
		cout << endl << "Loading world parameters... " << endl;
		
		vector<ID> id_au        = worldParams["id_au"];
		vector<string> name_au  = worldParams["name_au"];
		ID id_region            = worldParams["id_region"];
		string regionName       = worldParams["regionName"];
		
		vector<areaUnit> auvec = create_area_unit(id_au,
												  name_au,
												  id_region,
												  regionName);

		
		float unemployed_prop   = worldParams["unemployed_prop"];
		float pubT_prop         = worldParams["pubT_prop"];
		
		set_parameter(MP, "unemployed_prop", "double", worldParams);
		set_parameter(MP, "pubT_prop",       "double", worldParams);

		
		cout << "... world basic done ..." << endl;
		
		// Social place sizes distributions:
		
		vector<vector<uint> > size_hh_vec;
		vector<vector<uint> > size_wrk_vec;
		vector<vector<uint> > size_pubt_vec;
		vector<vector<uint> > size_school_vec;
		vector<vector<uint> > size_hosp_vec;
		vector<vector<uint> > size_other_vec;
		
		unsigned long n_au = id_au.size();
		
		List tmp_szdst;
		
		for(unsigned int i=0; i<n_au; i++)
		{
			List lst_hh		= worldParams["det_size_hh"];
			List lst_wkr	= worldParams["det_size_wrk"];
			List lst_pubt	= worldParams["det_size_pubt"];
			List lst_school = worldParams["det_size_school"];
			List lst_hosp	= worldParams["det_size_hosp"];
			List lst_other	= worldParams["det_size_other"];
			
			vector<uint>   hh_size		= lst_hh[i];
			vector<uint>   wrk_size		= lst_wkr[i];
			vector<uint>   pubt_size	= lst_pubt[i];
			vector<uint>   school_size	= lst_school[i];
			vector<uint>   hosp_size	= lst_hosp[i];
			vector<uint>   other_size	= lst_other[i];
			
			size_hh_vec.push_back(hh_size);
			size_wrk_vec.push_back(wrk_size);
			size_pubt_vec.push_back(pubt_size);
			size_school_vec.push_back(school_size);
			size_hosp_vec.push_back(hosp_size);
			size_other_vec.push_back(other_size);
		}
		
		cout << "... world deterministic SP sizes done ..." << endl;
		
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
		
		cout << "... all world parameters loaded." << endl;
		
		// === Schedules definition ===
		
		cout << "Creating schedules ..." << endl;
		// Define the time slices
		// must be same for all schedules and must sum up to 1.0
		vector<double> timeslice           = scheduleParams["timeslice"];
		vector<string> sched_names         = scheduleParams["sched_names"];
		vector<vector<string> > sched_desc = scheduleParams["sched_desc"];
		
		vector<schedule> sched = build_all_schedules(sched_desc, sched_names, timeslice);
		cout << "... schedules created." << endl;
		
		
		// =========================
		// === Call C++ function ===
		// =========================
		
		_RANDOM_GENERATOR.seed(rnd_seed);
		
		Simulator sim = run_detWorld(auvec,
									 size_hh_vec,
									 size_wrk_vec,
									 size_pubt_vec,
									 size_school_vec,
									 size_hosp_vec,
									 size_other_vec,
									 pr_age_hh,
									 unemployed_prop,
									 sched ,
									 MP,
									 start_time,
									 horizon,
									 i0,
									 interv_vec,
									 build_world_only);
		
		dcDataFrame world_sp  = sim.census_sp();
	
		Rcpp::List empty_list;
		
		bool debug_mode = MP.get_prm_bool("debug_mode");
		

		if(!build_world_only){
			
			// Retrieve all results from simulation:
			// populations:
			dcDataFrame world_final = export_world(sim.get_world(), light_output);
			
			// epidemic time series
			dcDataFrame ts = sim.timeseries();
			
			// Time series of census for every social places.
			//
			// Note: this is activated only in debug mode
			// because the exported objects are very large
			// and slow down everything .
			
			Rcpp::List census_sp;
			vector<string> tmp_names = {"time","id_sp","type","nS","nE"};
			
			if(debug_mode){
				census_sp.push_back(sim.get_ts_census_sp_time());
				census_sp.push_back(sim.get_ts_census_sp_id());
				census_sp.push_back(sim.get_ts_census_sp_type());
				census_sp.push_back(sim.get_ts_census_sp_nS());
				census_sp.push_back(sim.get_ts_census_sp_nE());
				census_sp.attr("names") = tmp_names;
			}
			else{
				census_sp = empty_list;
			}
			
			// Number of contacts tracker
			Rcpp::List track_n_contacts;
			track_n_contacts.push_back(sim.get_track_n_contacts_time());
			track_n_contacts.push_back(sim.get_track_n_contacts_uid());
			track_n_contacts.push_back(sim.get_track_n_contacts());
			tmp_names = {"time","uid","nContacts"};
			track_n_contacts.attr("names") = tmp_names;
			
			
			// Return R-formatted result:
			
			if(superlight_output){
				return List::create(Named("world_final")      = to_list(world_final,false),
									Named("ages")             = -999,
									Named("census_sp")		  = -999,
									Named("time_series")      = to_list(ts, false),
									Named("time_series_sp")   = census_sp,
									Named("track_n_contacts") = -999,
									Named("wiw_ages")         = -999,
									Named("contactAssort")    = -999,
									Named("sp_size_distrib")  = -999
									);

			}
			
			if(!superlight_output){
				return List::create(Named("world_final")      = to_list(world_final,false),
									Named("ages")             = census_ages(sim.get_world()),
									Named("census_sp")		  = to_list(world_sp,false),
									Named("time_series")      = to_list(ts, false),
									Named("time_series_sp")   = census_sp,
									Named("track_n_contacts") = track_n_contacts,
									Named("wiw_ages")         = sim.get_wiw_ages(),
									Named("contactAssort")    = sim.get_contactAssort(),
									Named("sp_size_distrib")  = sim.sp_size_distribution()
									);
			}
		}
		
		// If no simulation was requested,
		// return only information about the World created.
		if(build_world_only){
			
			// Return R-formatted result:
			return List::create(// DELETE WHEN SURE:Named("world")            = to_list_vector(W, false),
								Named("ages")             = census_ages(sim.get_world()),
								Named("census_sp")		  = to_list(world_sp,false),
								
								// kept for consistency with case
								// 'build_world_only=false':
								Named("population_final") = empty_list,
								Named("time_series")      = empty_list,
								Named("time_series_sp")   = empty_list,
								Named("track_n_contacts") = empty_list,
								Named("wiw_ages")         = empty_list,
								Named("contactAssort")    = empty_list
								);
		}
		
	}
	catch (...){
		::Rf_error(">>>> C++ exception (unknown reason) <<<<");
		return NULL;
	}
}


// [[Rcpp::export]]
List naf_run(List params,
             List simulParams,
             List intervParams,
             List worldParams,
             List scheduleParams,
			 unsigned int rnd_seed){
	
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
		set_parameter(MP, "contact_ratio_age_16_25",       "double", params);
		set_parameter(MP, "contact_ratio_age_25_40",       "double", params);
        set_parameter(MP, "contact_ratio_age_over_65",     "double", params);
        set_parameter(MP, "contact_ratio_sp_household",    "double", params);
        set_parameter(MP, "contact_ratio_sp_pubTransport", "double", params);
		
		set_parameter(MP, "contactAssort_lambda", "double", params);
        
        set_parameter(MP, "asymptom_infectiousness_ratio", "double", params);
        
        set_parameter(MP, "treat_doi_reduc", "double", params);
        set_parameter(MP, "treat_reduc_infect_mean", "double", params);
        
        set_parameter(MP, "vax_imm_hum_incr",          "double", params);
		set_parameter(MP, "vax_imm_cell_incr",          "double", params);
        set_parameter(MP, "vax_lag_full_efficacy", "double", params);
        
        set_parameter(MP, "proba_death_prm_1", "double", params);
        set_parameter(MP, "proba_death_prm_2", "double", params);
        set_parameter(MP, "proba_death_prm_3", "double", params);
        
		set_parameter(MP, "frailty_0",         "double", params);
		set_parameter(MP, "frailty_agepivot",  "double", params);
		set_parameter(MP, "frailty_slope1",    "double", params);
		set_parameter(MP, "frailty_slope2",    "double", params);
		set_parameter(MP, "frailty_sd",        "double", params);
		
		set_parameter(MP, "imm_hum_baseline",  "double", params);
		set_parameter(MP, "imm_hum_agezero",   "double", params);
		set_parameter(MP, "imm_hum_p",         "double", params);
		
		set_parameter(MP, "imm_cell_max",	"double", params);
		set_parameter(MP, "imm_cell_slope",	"double", params);
		set_parameter(MP, "imm_cell_pivot", "double", params);
        
        
        // ==== Simulator parameters ====
		
		bool build_world_only   = simulParams["build_world_only"];
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
			double	interv_efficacy     = tmp_interv["interv_efficacy"];
			double	interv_cvg_age_mult = tmp_interv["interv_cvg_age_mult"];
			
            intervention interv(interv_type,
                                interv_target,
                                interv_name,
                                interv_start,
                                interv_end,
                                interv_cvg_rate,
                                interv_cvg_max_prop,
								interv_efficacy,
								interv_cvg_age_mult);

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
		
		float unemployed_prop   = worldParams["unemployed_prop"];
		float pubT_prop         = worldParams["pubT_prop"];
		
		set_parameter(MP, "unemployed_prop", "double", worldParams);
		set_parameter(MP, "pubT_prop",       "double", worldParams);
		
        
        // Social place sizes distributions:
        
        vector<discrete_prob_dist<uint> > D_size_hh_vec;
        vector<discrete_prob_dist<uint> > D_size_wrk_vec;
        vector<discrete_prob_dist<uint> > D_size_pubt_vec;
        vector<discrete_prob_dist<uint> > D_size_school_vec;
        vector<discrete_prob_dist<uint> > D_size_hosp_vec;
        vector<discrete_prob_dist<uint> > D_size_other_vec;
        
        unsigned long n_au = id_au.size();
        
        List tmp_szdst;
        
        for(unsigned int i=0; i<n_au; i++)
		{
			
			vector<uint>   hh_size;
			vector<double> hh_size_proba;
			vector<uint>   wrk_size;
			vector<double> wrk_size_proba;
			vector<uint>   pubt_size;
			vector<double> pubt_size_proba;
			vector<uint>   school_size;
			vector<double> school_size_proba;
			vector<uint>   hosp_size;
			vector<double> hosp_size_proba;
			vector<uint>   other_size;
			vector<double> other_size_proba;
			
			read_discr_distr(worldParams,"hh_size",     hh_size, hh_size_proba, n_au, i);
			read_discr_distr(worldParams,"wrk_size",    wrk_size, wrk_size_proba, n_au, i);
			read_discr_distr(worldParams,"pubt_size",   pubt_size, pubt_size_proba, n_au, i);
			read_discr_distr(worldParams,"school_size", school_size, school_size_proba, n_au, i);
			read_discr_distr(worldParams,"hosp_size",   hosp_size, hosp_size_proba, n_au, i);
			read_discr_distr(worldParams,"other_size",  other_size, other_size_proba, n_au, i);
			
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
        
		Simulator sim = run_stochWorld(auvec,
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
									   unemployed_prop,
									   sched ,
									   MP,
									   start_time,
									   horizon,
									   i0,
									   interv_vec,
									   build_world_only);
		
		
		
		vector<dcDataFrame> W = export_dcDataFrame(sim.get_world());
		
		Rcpp::List empty_list;
		
		if(!build_world_only){
			
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
			
			
			
			// Return R-formatted result:
			return List::create(Named("population_final") = to_list(pop_final,false),
								// DELETE WHEN SURE:Named("world")            = to_list_vector(W, false),
								Named("time_series")      = to_list(ts, false),
								Named("time_series_sp")   = census_sp,
								Named("track_n_contacts") = track_n_contacts,
								Named("wiw_ages")         = sim.get_wiw_ages(),
								Named("contactAssort")    = sim.get_contactAssort()
								);
		}
		
		// If no simulation was requested,
		// return only information about the World created.
		if(build_world_only){
			
			dcDataFrame world_sp = sim.census_sp();
			
			// Return R-formatted result:
			return List::create(// DELETE WHEN SURE:Named("world")            = to_list_vector(W, false),
								Named("ages")             = census_ages(sim.get_world()),
								Named("census_sp")		  = to_list(world_sp,false),
								
								// kept for consistency with case
								// 'build_world_only=false':
								Named("population_final") = empty_list,
								Named("time_series")      = empty_list,
								Named("time_series_sp")   = empty_list,
								Named("track_n_contacts") = empty_list,
								Named("wiw_ages")         = empty_list,
								Named("contactAssort")    = empty_list
								);
		}
		
    }
    catch (...){
        ::Rf_error(">>>> C++ exception (unknown reason) <<<<");
        return NULL;
    }
}


