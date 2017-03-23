
######    FUNCTION UTILITIES TO RUN SIMULATIONS



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




earliest.interv.start <- function(L){
	min.start <- 9999999
	for(i in seq_along(L)){
		a <- L[[i]]
		for(j in seq_along(a)) {
			tmp <- a[[j]][['interv_start']]
			if(tmp <min.start) min.start <- tmp
		}
	}
	return(min(-1, min.start))
}



# Run the simulation for a given 
# set of scenario parameters
run.simul <- function(scen.id, 
					  dir.save.rdata = './', 
					  baseonly = FALSE, 
					  force.light.output = FALSE) {
	
	t0        <- as.numeric(Sys.time())
	n.MC      <- simul.prm[['mc']]
	n.cpu     <- simul.prm[['cpu']]
	seeds     <- 1:n.MC
	
	# optimize simulation start time based 
	# on interventions:
	L <- list(interv.prm.0, interv.prm)
	old.start <- simul.prm[['start_time']]
	simul.prm[['start_time']] <- earliest.interv.start(L) -2 
	print(paste('Start time optimized from ',old.start, ' to ', simul.prm[['start_time']] ))
	
	if(force.light.output) simul.prm[['light_output']] <- TRUE
	
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
	res.list <- NA
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
	
	if(scen.id > 0){
	filename <- paste0(dir.save.rdata,'mc-simul-',scen.id,'.RData')

	}
	if(scen.id<= 0) {
		filename <- paste0(dir.save.rdata,'mc-simul.RData')
		# save.image(file= filename, compress = FALSE,envir = environment())
	}
	
	save(list=c('res.list.0', 'res.list', 
				'simul.prm',  'world.prm'), 
		 file= filename, compress = FALSE)
	
	print('... RData file saved.')
	t1 <- as.numeric(Sys.time())
	print(paste("Scenario",scen.id,": Simulation computing time is",round((t1-t0)/60,1),"min"))
}



# Run the simulation. Used to fit R
run.simul.fit.R <- function(scen.id, dir.save.rdata, prm, simul.prm, interv.prm.0, world.prm, sched.prm, stoch_build_world) {
	
    t0        <- as.numeric(Sys.time())
    n.MC      <- simul.prm[['mc']]
    n.cpu     <- simul.prm[['cpu']]
    seeds     <- 1:n.MC
    
    # optimize simulation start time based 
    # on interventions:
    L <- list(interv.prm.0, interv.prm)
    old.start <- simul.prm[['start_time']]
    simul.prm[['start_time']] <-  -2 
    print(paste('Start time forced to',simul.prm[['start_time']]))
    
	sfInit(parallel = (n.cpu>1), cpu = n.cpu)
	sfLibrary(naf,lib.loc = R.library.dir) 
	
	# Baseline scenario only:
	res.list.0 <- sfSapply(seeds, run.snow.wrap,
						   prm        = prm, 
						   simul.prm  = simul.prm, 
						   interv.prm = interv.prm.0, 
						   world.prm  = world.prm, 
						   sched.prm  = sched.prm,
						   stoch_build_world = stoch_build_world,
						   simplify   = FALSE)
	
	sfStop()
	filename <- paste0(dir.save.rdata,'mc-simul-',scen.id,'.RData')
	save(list=c('res.list.0', 'simul.prm',  'world.prm'), 
	     file= filename, compress = FALSE)
	
	print('... RData file saved.')
	t1 <- as.numeric(Sys.time())
	print(paste("Scenario",scen.id,": Simulation computing time is",round((t1-t0)/60,1),"min"))
}



