//
//  build_world.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-08-02.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "build_world.h"


vector<individual> create_individuals(uint n, uint first_id_indiv){
    /// Build several (empty) individuals
    
    vector<individual> x(n);
    for (int i=0; i<n; i++){
        individual tmp(first_id_indiv + i);
        x[i] = tmp;
    }
    return x;
}


void keep_indiv_with_household(vector<individual>& x){
    /// Keep individuals in this vector that are linked to a household.
    /// Individuals who are _not_ linked, are deleted.
    
    cout << endl<<" Removing individuals w/o households..."; fflush(stdout);
    
    unsigned long initial_size = x.size();
    uint cnt = 0;
    for (uint i=0; i<x.size(); i++) {
        if (x[i].get_id_sp_household() == __UNDEFINED_ID) {
            x.erase(x.begin()+i);
            i--;
            cnt++;
        }
    }
    cout << " done: "<<cnt<< " individuals removed out of a total of ";
    cout << initial_size<<" ("<< (int)((double)(cnt)/initial_size*100) <<"%)."<<endl;
}



void keep_indiv_with_relevant_links(vector<individual>& x){
    /// Keep individuals in this vector that are linked to a scheduled SP type.
    /// Individuals who are _not_ linked, are deleted.

    cout << endl<<" Removing individuals w/o link to scheduled SP..."; fflush(stdout);
    
    unsigned long initial_size = x.size();
    uint cnt = 0;
    for (uint i=0; i<x.size(); i++) {
        if (x[i].get_id_sp_household() == __UNDEFINED_ID) {
            x.erase(x.begin()+i);
            i--;
            cnt++;
        }
        else if(x[i].get_id_sp_hospital() == __UNDEFINED_ID){
            x.erase(x.begin()+i);
            i--;
            cnt++;
        }
        else if(x[i].get_id_sp_other() == __UNDEFINED_ID){
            x.erase(x.begin()+i);
            i--;
            cnt++;
        }
        else if(x[i].get_id_sp_pubTransp() == __UNDEFINED_ID){
            x.erase(x.begin()+i);
            i--;
            cnt++;
        }
        else if(x[i].get_id_sp_school() == __UNDEFINED_ID){
            x.erase(x.begin()+i);
            i--;
            cnt++;
        }
        else if(x[i].get_id_sp_workplace() == __UNDEFINED_ID){
            x.erase(x.begin()+i);
            i--;
            cnt++;
        }
    }
    cout <<" done: "<<cnt<< " individuals removed out of a total of "<<initial_size;
    cout <<" ("<< (int)((double)(cnt)/initial_size*100) <<"%)."<<endl;
}



vector<areaUnit> create_area_unit(const vector<ID>& id_au,
                                  const vector<string>& name_au,
                                  ID id_region,
                                  string regionName){
    vector<areaUnit> v;
    unsigned long n = id_au.size();
    
    stopif(n != name_au.size(), "vector size do not match.");
    
    for (uint i=0; i<n; i++) {
        areaUnit tmp(id_au[i], name_au[i], id_region, regionName);
        v.push_back(tmp);
    }
    return v;
}


vector<socialPlace> create_socialPlaces_size(SPtype sp_type,
                                             uint num_sp,
                                             uint first_id_sp,
                                             uint first_id_indiv,
                                             discrete_prob_dist<uint> size_distrib,
                                             areaUnit AU,
                                             vector<individual>& indiv,
                                             float age_min, float age_max){
    /// Build 'num_sp' social of 'sp_type' according
    /// to a pre-specified size distribution.
    /// Individuals passed in parameters of this function
    /// are linked to each social place.
    /// Note: There is a constraint on individual's age to be linked.
    /// * * WARNING * *
    /// vector of individuals must be large enough!
    
//    cout << endl <<  "Creating social places of type " + SPtype2string(sp_type) + " ...";
//    fflush(stdout);
    
    vector<socialPlace> x;
    unsigned long n_indiv = indiv.size();
    
    uint seed = 1234;
    
    // Draw the size for each social place:
    vector<uint> size_sp = size_distrib.sample(num_sp, seed);
    
    uint cnt = 0; // <--- counts over the _whole_ vector of individuals supplied and make sure individual id do not overlap across AUs.
    uint sp_not_linked = 0;
    
    for (ID k=0; k < num_sp; k++) {
        
        // Create social place (empty shell):
        socialPlace tmp(AU, first_id_sp + k, sp_type);
        
        bool record = false;
        
        // Link individuals with the kth social place:
        for (uint linked = 0;
             linked < size_sp[k] &&
             cnt < n_indiv; )
        {
            float age = indiv[cnt].get_age(); // <-- no 'first_id_indiv' for the index of vector ('first_id_indiv' is just for ID)
            
            // DEBUG
            //cout << SPtype2string(sp_type) << ": SP_"<<first_id_sp + k <<"  ; first_id_sp = "<<first_id_sp << " k = "<<k << endl;
            
            if(age_min <= age && age<age_max)
            {
                // Link both individual and social place
                // ('set_id_sp' does both):
                indiv[cnt].set_id_sp(sp_type, tmp);
                record = true;
                linked++;
            }
            cnt++;
        }
        // Record the links
        if(record)
            x.push_back(tmp);

        if(cnt == n_indiv){
            // end of indiv vector reached
            sp_not_linked++;
        }
        
    }
    if(sp_not_linked>0){
        cout << endl;
        cout << "* * WARNING * * \n ";
        cout << sp_not_linked;
        cout << " out of the "<< num_sp <<" ";
        cout << SPtype2string(sp_type) << " social places requested";
        cout <<" were not created/completed because of lack of individuals.";
        cout << endl;
    }
    
    cout << "Social places of type " + SPtype2string(sp_type) + " completed." <<endl;
    return x;
}



vector<socialPlace> create_other_socialPlaces(uint num_sp,
                                              uint first_id_sp,
                                              discrete_prob_dist<uint> size_distrib,
                                              areaUnit AU){
    /// Create social places of type 'other'.
    /// This type of social place does _not_ have linked individual.
    /// (they will come here randomly)
    
    vector<socialPlace> x;
    uint seed = 1234;
    
    // Draw the size for each social place:
    vector<uint> size_sp = size_distrib.sample(num_sp, seed);
    
    for (ID k=0; k < num_sp; k++) {
        
        // Create social place (empty shell):
        socialPlace tmp(AU, first_id_sp + k, SP_other);
        x.push_back(tmp);
    }
    return x;
}




void assign_age_in_households(vector<socialPlace>& hh,
                              vector<individual>& indiv,
                              vector<vector<discrete_prob_dist<uint> > > age_distrib){
    /// Assign the age of individuals based on
    /// the size of the household they are in.
    
    uint seed = 0;
    
    unsigned long n_ad = age_distrib.size();
    
    for (uint k=0; k<hh.size(); k++)
    {
        uint hh_size = hh[k].n_linked_indiv();
        uint idx     = hh_size -1; // <-- warning '-1'
        
        stopif(idx >= n_ad,
               "Age distribution within households needed (but not defined) for sizes larger than " + to_string(n_ad)+".");
        
        for (uint i=0; i<hh_size; i++) {            
            // Retrieve the appropriate age distribution
            // for the size of this social place:
            discrete_prob_dist<uint> pr = age_distrib[idx][i];
            uint age = pr.sample(1, seed++)[0];
            
            // Assign the age to the linked individual:

            // 1) Retrieve the ID of the ith individual in this household
            ID id_tmp = hh[k].get_linked_indiv_id(i);
            
            // 2) Retrieve its position in the 'indiv' vector
            uint pos  = 0;
            while (indiv[pos].get_id() != id_tmp && pos<indiv.size()) {
                pos++;
            }
            stopif(pos==indiv.size(),"Individual ID "+ to_string(id_tmp) +" not found.");

            // 3) assign age
            indiv[pos].set_age(age);
        }
    }
}


void populate_households(vector<socialPlace>& sp_hh,
                         vector<individual>& indiv){
    /// Populate (i.e. fill vector '_indiv') with
    /// individuals that are linked to this social place.
    
    for (uint k=0; k<sp_hh.size(); k++) {
        uint hh_size = sp_hh[k].n_linked_indiv();
        for (uint i=0; i<hh_size; i++) {
            // Assign the age to the linked individual:
            ID id_tmp = sp_hh[k].get_linked_indiv_id(i);
            individual tmp = get_indiv_with_ID(id_tmp, indiv);
            sp_hh[k].add_indiv(tmp);
        }
    }
}



vector<socialPlace> build_world(vector<areaUnit> AU,
                                vector<discrete_prob_dist<uint> > D_size_hh,     // Households sizes
                                vector< vector<discrete_prob_dist<uint> > > pr_age_hh,  // Age distribution inside households
                                
                                vector<discrete_prob_dist<uint> > D_size_wrk,
                                vector<discrete_prob_dist<uint> > D_size_pubt,
                                vector<discrete_prob_dist<uint> > D_size_school,
                                vector<discrete_prob_dist<uint> > D_size_hosp,
                                vector<discrete_prob_dist<uint> > D_size_other,
                            
                                vector<uint> n_hh,
                                vector<uint> n_wrk,
                                vector<uint> n_pubt,
                                vector<uint> n_school,
                                vector<uint> n_hosp,
                                vector<uint> n_other){
    /// Build world before simulation runs
    
    vector< vector<socialPlace> > y;
    
    // Variables that make sure id (sp & indiv)
    // do not overlap across Area Units!
    uint first_id_sp       = 0; // make sure no overlap for id_sp _within_ an Area Unit
    uint first_id_indiv_au = 0; // make sure no overlap for id_indiv _across_ Area Units
    
    for (uint a=0; a<AU.size(); a++)
    {
        // These individuals are the raw material
        // for creating the world.
        
        // First, make sure there are enough individuals created,
        // by taking the maxing value of the support of the
        // distribution of households sizes:
        vector<uint> ds          = D_size_hh[a].get_value();
        uint nmax                = *std::max_element(ds.begin(), ds.end());
        vector<individual> indiv = create_individuals(n_hh[a]*nmax, first_id_indiv_au);
        
        cout << indiv.size() << " individuals provided to build area unit ID_"<< a << endl;

        // Create households and link existing individuals to them.
        // The size of the households is driven by the
        // distribution 'D_hh_size'.
        vector<socialPlace> sp_hh   = create_socialPlaces_size(SP_household,
                                                               n_hh[a],
                                                               first_id_sp,
                                                               first_id_indiv_au,
                                                               D_size_hh[a],
                                                               AU[a],
                                                               indiv);
        // Assign individuals' age based on the size of
        // the household they are linked to.
        assign_age_in_households(sp_hh, indiv, pr_age_hh);
        
        // Get rid of excess individuals:
        // (because not all of them could
        //  be allocateda scheduled SP)
        keep_indiv_with_household(indiv);
        
        // Make sure ID of social places do not
        // overlap within _and_ across area units:
        first_id_sp += (uint) sp_hh.size();
        
        
        // Link social places (other than households)
        // to individuals.
        // Constraints to be linked to a given social place
        // is based on
        // - the distribution of sizes for this social place
        // - individual's age
        
        // only indiv older than 18 yrs old
        // can be linked to workplaces:
        float age_min_wrk = 17.9;
        float age_max_wrk = _AGE_MAX; // <-- TO DO: change that!
        
        vector<socialPlace> sp_wrk = create_socialPlaces_size(SP_workplace,
                                                              n_wrk[a],
                                                              first_id_sp,
                                                              first_id_indiv_au,
                                                              D_size_wrk[a],
                                                              AU[a],
                                                              indiv,
                                                              age_min_wrk,
                                                              age_max_wrk);
        
        first_id_sp += (uint) sp_wrk.size()  ;
        
        float age_min_pubt = 0.0;
        float age_max_pubt = _AGE_MAX;
        vector<socialPlace> sp_pubt = create_socialPlaces_size(SP_pubTransp,
                                                               n_pubt[a],
                                                               first_id_sp,
                                                               first_id_indiv_au,
                                                               D_size_pubt[a],
                                                               AU[a],
                                                               indiv,
                                                               age_min_pubt,
                                                               age_max_pubt);
        first_id_sp += (uint) sp_pubt.size() ;
        
        // Schools are for children only.
        // In particular, teachers and university students
        // are not modelled with social place school.
        // - For teachers, they are allocated to a workplace
        // (although they do not mix as much as children b/w
        //  them, that would be better to include them in school --> TO DO)
        // - University students are linked to a standard workplace, this
        //  is justified by the fact they tend to mix more like adults at that age.
        
        float age_min_school = 0.0;
        float age_max_school = 18.1;
        vector<socialPlace> sp_school = create_socialPlaces_size(SP_school,
                                                                 n_school[a],
                                                                 first_id_sp,
                                                                 first_id_indiv_au,
                                                                 D_size_school[a],
                                                                 AU[a],
                                                                 indiv,
                                                                 age_min_school,
                                                                 age_max_school);
        // No age constraint for linking individuals
        // to hospital and 'other' type of social places.
       
        
        first_id_sp += (uint) sp_school.size();
        
        vector<socialPlace> sp_hosp = create_socialPlaces_size(SP_hospital,
                                                               n_hosp[a],
                                                               first_id_sp,
                                                               first_id_indiv_au,
                                                               D_size_hosp[a],
                                                               AU[a],
                                                               indiv);
        
        first_id_sp += (uint) sp_hosp.size() ;
        
        vector<socialPlace> sp_other = create_socialPlaces_size(SP_other,
                                                                n_other[a],
                                                                first_id_sp,
                                                                first_id_indiv_au,
                                                                D_size_other[a],
                                                                AU[a],
                                                                indiv);
        
        // Makes sure sp_ids are incremented
        // for the next AU iteration:
        first_id_sp += (uint) sp_other.size() ;
        
        // 'Physically' move individuals
        // in the household they are linked to:
        populate_households(sp_hh, indiv);
        
        // Put all social places of
        // the area unit into one container:
        vector<vector<socialPlace> > tmp;
        tmp.push_back(sp_hh);
        tmp.push_back(sp_wrk);
        tmp.push_back(sp_pubt);
        tmp.push_back(sp_school);
        tmp.push_back(sp_hosp);
        tmp.push_back(sp_other);
        
        // Transform the container into
        // one tall vector:
        vector<socialPlace> tmp2 = melt(tmp);
        check_sp_integrity(tmp2);

        // Increase the starting value of IDs for
        // individuals of the next area unit:
        first_id_indiv_au += indiv.size();
        
        // store the tall containers for each area unit:
        y.push_back(tmp2);
    }
    // Transform the container of all social places in
    // all area units into one tall vector:
    vector<socialPlace> x = melt(y);
    check_sp_integrity(x);
    
    return x;
}


void check_sp_integrity(const vector<socialPlace>& x){
    
    vector<uint> id_sp(x.size());
    for(uint i=0;i<x.size();i++) id_sp[i]=x[i].get_id_sp();
    vector<uint> tmp = sort_remove_duplicate(id_sp, true);
    stopif(tmp.size()<id_sp.size(), "SP IDs not uniques");

    vector<uint> id_indiv;
    for(uint i=0;i<x.size();i++) {
        for(uint j=0;j<x[i].get_indiv().size();j++){
            id_indiv.push_back(x[i].get_indiv(j).get_id());
        }
    }
    tmp.clear();
    tmp = sort_remove_duplicate(id_indiv, true);
    
    // DEBUG
    if(tmp.size()<id_indiv.size()){
        displayVector(id_indiv);
        displayVector(tmp);
    }
    stopif(tmp.size()<id_indiv.size(), "individual IDs not uniques");
}



void assign_schedules(vector<socialPlace> & W,
                      const vector<schedule>& sched,
                      float prop_unemployed){
    /// Assign schedule for all individual in the world.
    /// For each individual, schedule chosen based on
    /// age, employment status, etc.
    
    
    uint s_student      = pos_schedname("student", sched);
    uint s_worker_sed   = pos_schedname("worker_sed", sched);
    uint s_unemployed   = pos_schedname("unemployed", sched);
    //uint s_worker_trav  = pos_schedname("worker_trav", sched);
    
    for (uint k=0; k<W.size(); k++)
    {
        unsigned long nk = W[k].get_size();
        
        for (uint i=0; i< nk; i++)
        {
            float age = W[k].get_indiv(i).get_age();

            if (0.0 <= age && age < 18.0) {
                W[k].set_schedule_indiv(i, sched[s_student]);
            }
            else if (age >= 18.0){
                
                // Unemployed
                std::uniform_real_distribution<> unif01(0.0, 1.0);
                double u = unif01(_RANDOM_GENERATOR);
                if (u < prop_unemployed){
                    W[k].set_schedule_indiv(i, sched[s_unemployed]);
                }
                else{
                    W[k].set_schedule_indiv(i, sched[s_worker_sed]);
                }
            }
            bool check = (W[k].get_indiv(i).get_schedule().get_sp_type().size()==0);
            if(check) W[k].get_indiv(i).displayInfo();
            stopif(check, "Schedule not assigned for indiv_"+ to_string(i)+ " (currently in SP_"+to_string(k)+").");
        }
    }
}















