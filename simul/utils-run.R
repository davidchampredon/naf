
##################################################################
######
######    FUNCTION UTILITIES TO RUN SIMULATIONS
######
##################################################################



# Snowfall wrapper
run.snow.wrap <- function(seedMC,
						  prm, 
						  simul.prm, 
						  interv.prm, 
						  world.prm, 
						  sched.prm,
						  stoch_build_world){
	
	if(stoch_build_world){
		res <- naf_run(prm, 
					   simul.prm, 
					   interv.prm, 
					   world.prm, 
					   sched.prm,
					   seedMC)
	}
	else{
		res <- naf_run_det(prm, 
						   simul.prm, 
						   interv.prm, 
						   world.prm, 
						   sched.prm,
						   seedMC)
	}
	return(res)
}

# Run the simulation for a given 
# set of scenario parameters
run.simul <- function(scen.id, dir.save.rdata = './') {
	
	t0 <- as.numeric(Sys.time())
	baseonly  <- FALSE # <-- force scneario comparison   ((simul.prm[['baseline_only']]))
	n.MC      <- simul.prm[['mc']]
	n.cpu     <- simul.prm[['cpu']]
	seeds     <- 1:n.MC
	
	sfInit(parallel = (n.cpu>1), cpu = n.cpu)
	sfLibrary(naf,lib.loc = R.library.dir) 
	
	# Baseline scenario:
	res.list.0 <- sfSapply(seeds, run.snow.wrap,
						   prm        = prm, 
						   simul.prm  = simul.prm, 
						   interv.prm = interv.prm.0, 
						   world.prm  = world.prm, 
						   sched.prm  = sched.prm,
						   stoch_build_world = stoch_build_world,
						   simplify   = FALSE)
	
	# Interventions:
	if(!baseonly){
		print('Starting simulations with interventions...')
		res.list <- sfSapply(seeds, run.snow.wrap,
							 prm        = prm, 
							 simul.prm  = simul.prm, 
							 interv.prm = interv.prm, 
							 world.prm  = world.prm, 
							 sched.prm  = sched.prm,
							 stoch_build_world = stoch_build_world,
							 simplify   = FALSE)
		print('... simulations with interventions done.')
	}
	sfStop()
	
	print('Saving RData file...')
	save(list=c('res.list.0','res.list'), 
		 file= paste0(dir.save.rdata,'mc-simul-',scen.id,'.RData'),
		 compress = FALSE)
	print('... RData file saved.')
	t1 <- as.numeric(Sys.time())
	print(paste("Scenario",scen.id,": Simulation computing time is",round((t1-t0)/60,1),"min"))
}
