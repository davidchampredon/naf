##################################################################
######
######    RUN SIMULATION ON ONE OR TWO SCENARIOS
######
######    
######
######
##################################################################

t0 <- as.numeric(Sys.time())

R.library.dir    <- '../Rlibrary/lib'
param.model.dir  <- '../param-model/'
data.dir         <- '../data/'

library(ggplot2)
library(plyr)
library(naf,lib.loc = R.library.dir)
library(snowfall)

source('analysis_tools.R')
source('build-world-det-sizegen.R')
source('utils-run.R')
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

save.plot.to.file <- TRUE

prm        <- list()
simul.prm  <- list()
interv.prm <- list()
world.prm  <- list()
sched.prm  <- list()


### ==== Model Parameters ====

prm[['debug_mode']] <- FALSE

prm       <- read.prm(filename = paste0(param.model.dir,fname.prm.epi), 
					  prm = prm)
simul.prm <- read.prm(filename = paste0(param.model.dir,fname.prm.simul), 
					  prm = simul.prm)

stoch_build_world <- simul.prm[['build_world_stoch']]

###  ==== World parameters ====

world.prm <- load.world.prm(filename = paste0(param.model.dir,fname.prm.au), 
							path.prm = param.model.dir)

world.prm[['id_region']]  <- 0
world.prm[['regionName']] <- "Canada"
world.prm[['unemployed_prop']] <- 0.10

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

# schedule time slices:

sched.prm[['timeslice']] <- c(1.0/24, 4.0/24, 4.0/24, 
							  1.0/24, 2.0/24, 12.0/24)

sched.prm[['sched_desc']]  <- read.schedules(paste0(param.model.dir,fname.schedules), return.names = FALSE)
sched.prm[['sched_names']] <- read.schedules(paste0(param.model.dir,fname.schedules), return.names = TRUE)

###  ==== Intervention parameters ====

# Baseline, no vax intervention:
interv.prm.0 <- load.interv.prm(paste0(param.model.dir,fname.prm.interv.0))
# Vax intervention:
interv.prm   <- load.interv.prm(paste0(param.model.dir,fname.prm.interv))


### === Snowfall wrapper ===

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


### ==== Run Simulation ====

baseonly  <- simul.prm[['baseline_only']]

run.simul(scen.id = -1, 
		  dir.save.rdata = './',
		  baseonly = baseonly) 

