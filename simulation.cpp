//
//  simulation.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-06.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "simulation.h"




void move_indiv(const SPtype sptype, world& spvec, double proba){
	
	/// Move individuals across social places
	
	for (int k=0; k<spvec.size(); k++)
	{
		vector<unsigned int> pos2move; // idx position to move (must be done at the end of the loop, bc vector keeps on changing size!)
		for (int i=0; i<spvec[k].get_size(); i++)
		{
			// Draw the chance move will actually happen:
			double u = (double)(rand()) / RAND_MAX;
			if ( u < proba )
			{
				individual tmp = spvec[k].get_indiv()[i];
				
				ID id_dest = __UNDEFINED_ID;
				
				if(sptype == SP_household)	id_dest = tmp.get_id_sp_household();
				if(sptype == SP_workplace)	id_dest = tmp.get_id_sp_workplace();
				if(sptype == SP_school)		id_dest = tmp.get_id_sp_school();
				if(sptype == SP_other)		id_dest = tmp.get_id_sp_other();
				if(sptype == SP_hospital)	id_dest = tmp.get_id_sp_hospital();
				if(sptype == SP_pubTransp)	id_dest = tmp.get_id_sp_pubTransp();
				
				if(id_dest != __UNDEFINED_ID &&
				   id_dest != k)
				{
					// add indiv to destination
					spvec[id_dest].add_indiv(tmp);
					// record its position for future deletion
					pos2move.push_back(i);
				}
				// remove from this SP the individuals that moved:
				if(pos2move.size()>0) spvec[k].remove_indiv(pos2move);
			}
		}
	}
}


unsigned int transmission(socialPlace& sp, double contact_rate, double dt){
	
	/// Performs transmission within a social place. Returns incidence
	
	
	unsigned int n = sp.get_prevalence();
	unsigned int inc = 0;
	
	if (n > 0){
		
		// Calculate total number of contacts:
		unsigned int nContacts = (unsigned int)(contact_rate * n * dt);

		cout << " DEBUG TRANSMISSION: nContacts: " << nContacts << endl;
		
		// number of susceptible in this social place:
		//unsigned int nS = sp.get_size() - n;
		
		// Susceptible candidates for transmission:
		vector<unsigned int> susc = sp.pick_rnd_susceptibles(nContacts);
		
		
		cout << " DEBUG TRANSMISSION: #susc: " << susc.size() << endl;
		
		// Calculate transmission based on susceptible features:

		for (unsigned int i=0; i<susc.size(); i++) {

			// Transmission Probability
			
			double imm = sp.get_indiv()[susc[i]].get_immunity();
			double fra = sp.get_indiv()[susc[i]].get_frailty();
			double pt = (1-imm) * fra ;
			
			cout << " DEBUG TRANSMISSION: pt: " << pt << endl;

			
			// Draw event:
			double u = (double)(rand())/RAND_MAX;
			if( u < pt) {
				sp.get_indiv()[susc[i]].acquireDisease();
				sp.increase_prevalence();
				inc ++;
			}
		}
	}
	return inc;
}