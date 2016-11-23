##################################################################
######
######    UTILITIES SCRIPT TO LOAD SIMULATION PARAMETERS
######
######    - Source this file to load all parameters for one scenario.
######    - Several parameters lists will be created.
######    - Then overwrite parameters values if several scenatios
######    are requested.
######
##################################################################

library(ggplot2)
library(plyr)
library(snowfall)

try({
	R.library.dir    <- '../Rlibrary/lib'
	param.model.dir  <- '../param-model/'
	data.dir         <- '../data/'
	
	library(naf,lib.loc = R.library.dir)
	
	source('analysis_tools.R')
	source('build-world-det-sizegen.R')
	source(paste0(param.model.dir,'read-prm.R'))
	source(paste0(param.model.dir,'read-world-prm.R'))
	source(paste0(param.model.dir,'read-interv.R'))
	source(paste0(param.model.dir,'read-sched.R'))
	
	
	# Parameter file names
	fname.prm.epi      <- 'prm-epi.csv'
	fname.prm.simul    <- 'prm-simul.csv'
	fname.prm.au       <- 'prm-au-ontario.csv'
	fname.prm.interv.0 <- 'prm-interv-0.csv'
	fname.prm.interv   <- 'prm-interv.csv'
	fname.schedules    <- 'schedules.csv'
	
	prm        <- list()
	simul.prm  <- list()
	interv.prm <- list()
	world.prm  <- list()
	sched.prm  <- list()
	
	
	### ==== Model Parameters ====
	
	prm[['debug_mode']] <- FALSE
	
	prm       <- read.prm(filename = paste0(param.model.dir, fname.prm.epi), 
						  prm = prm)
	simul.prm <- read.prm(filename = paste0(param.model.dir, fname.prm.simul), 
						  prm = simul.prm)
	
	stoch_build_world <- simul.prm[['build_world_stoch']]
	
	
	###  ==== World parameters ====
	
	world.prm <- load.world.prm(filename = paste0(param.model.dir, fname.prm.au), 
								path.prm = param.model.dir)
	
	world.prm[['id_region']]  <- 0
	world.prm[['regionName']] <- "Canada"
	world.prm[['unemployed_prop']] <- 0.10
	
	# rescale world size:
	sf           <- as.numeric(simul.prm[['scale_factor']])
	world.prm    <- scale.world(1/sf, world.prm)
	print(paste('World size reduced by',sf))
	
	# age distributions, conditional on 
	# household composition:
	base.fname  <- paste0(data.dir,'households/hh_size_ad')
	world.prm   <- c(world.prm,
					 load.age.hh(base.fname, 
					 			max.hh.size = max(world.prm[['max_hh_size']])) )
	
	world.prm <- gen_world_ontario(world.prm)
	
	
	#  ==== Schedules ====
	
	sched.prm[['timeslice']]   <- c(1.0/24, 4.0/24, 4.0/24, 
									1.0/24, 2.0/24, 12.0/24)
	sched.prm[['sched_desc']]  <- read.schedules(paste0(param.model.dir,fname.schedules), return.names = FALSE)
	sched.prm[['sched_names']] <- read.schedules(paste0(param.model.dir,fname.schedules), return.names = TRUE)
	
	
	###  ==== Intervention parameters ====
	
	# Baseline, no vax intervention:
	interv.prm.0 <- load.interv.prm(paste0(param.model.dir,fname.prm.interv.0))
	# Vax intervention:
	interv.prm   <- load.interv.prm(paste0(param.model.dir,fname.prm.interv))
	
},
silent = TRUE)

overwrite.selected.param <- function(filename, scen.id) {
	
	x <- read.csv(filename)
	i <- which(x$scenario_id == scen.id)
	y <- x[i,]
	
	# WARNING:
	# interv.prm[[1]] means the FIRST intervention is overwritten.
	
	interv.prm[[1]][['interv_target']]   <<- as.character(y$interv_target)
	interv.prm[[1]][['interv_start']]    <<- y$interv_start
	interv.prm[[1]][['interv_cvg_rate']] <<- y$interv_cvg_rate
	
	prm[['contact_rate_mean']]   <<- y$contact_rate_mean
	prm[['contact_rate_stddev']] <<- y$contact_rate_stddev
	
	prm[['vax_imm_hum_incr']]    <<- y$vax_imm_hum_incr
	prm[['vax_imm_cell_incr']]   <<- y$vax_imm_cell_incr
	
	prm[['imm_hum_baseline']]    <<- y$imm_hum_baseline
}



