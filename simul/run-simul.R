##################################################################
######
######    TEST FOR 'NAF' LIBRARY
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

source('analysis_tools.R')
source(paste0(param.model.dir,'read-prm.R'))
source(paste0(param.model.dir,'read-world-prm.R'))
source(paste0(param.model.dir,'read-interv.R'))

# Parameter file names
fname.prm.epi    <- 'prm-epi.csv'
fname.prm.simul  <- 'prm-simul.csv'
fname.prm.au     <- 'prm-au-ontario.csv'
fname.prm.interv <- 'prm-interv.csv'

do.plot           <- TRUE 
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

###  ==== World parameters ====

world.prm <- load.world.prm(filename = paste0(param.model.dir,fname.prm.au), 
							path.prm = param.model.dir)

world.prm[['id_region']]  <- 0
world.prm[['regionName']] <- "Canada"

sf           <- as.numeric(simul.prm[['scale_factor']])
world.prm    <- scale.world(1/sf, world.prm)
print(paste('World size reduced by',sf))

# age distributions, conditional on 
# household composition:
base.fname  <- paste0(data.dir,'households/hh_size_ad')
world.prm   <- c(world.prm,
				 load.age.hh(base.fname, 
				 			 max(world.prm[['max_hh_size']])) )

# schedule time slices:
sched.prm[['timeslice']] <- c(1.0/24, 4.0/24, 4.0/24, 
							  1.0/24, 2.0/24, 12.0/24)


###  ==== Intervention parameters ====

interv.prm <- load.interv.prm(paste0(param.model.dir,fname.prm.interv))


### ==== Run Simulation ====

res.try <- try(
	res <- naf_run(prm, 
				   simul.prm, 
				   interv.prm, 
				   world.prm, 
				   sched.prm)
)

if(class(res.try)=='try-error') stop()

t1 <- as.numeric(Sys.time())


# ==== Process results ====

message("Processing results...")

ts <- as.data.frame(res[['time_series']])

# JUST ONE SOCIAL PLACE--->   pop <- as.data.frame(res[['population_final']])
world0 <- res[['world']]
z      <- lapply(world0, as.data.frame)
pop    <- do.call('rbind',z)

ws <- ddply(pop, c('id_au','sp_type'), summarize, 
			n_sp    = length(id_sp),
			n_indiv = length(id_indiv))
message("Processing done.")


### ==== PLOTS ==== 

if(do.plot){
	message("Plotting results...")
	
	# Population:
	if (save.plot.to.file) pdf('plot_pop.pdf', width = 30,height = 18)
	try( plot.population(pop),  silent = T)
	try( plot.n.contacts(res$track_n_contacts),  silent = T)
	try( plot.age.contact.matrix(res),  silent = T)
	try( plot.sp.sz.distrib(pop,world.prm) , silent = T)
	try( plot.share.same.hh(pop), silent = F)
	if (save.plot.to.file) dev.off()
	
	# Time series:
	if (save.plot.to.file) pdf('plot_ts.pdf', width = 25,height = 15)
	
	try( plot.epi.timeseries(ts),  silent = T)
	try( grid.arrange(plot.ts.sp(res$time_series_sp),
					  plot.ts.sp(res$time_series_sp, facets = T)), 
		 silent = T)
	
	if (save.plot.to.file) dev.off()
	message("Plotting done.")
}
	
# End
t2 <- as.numeric(Sys.time())
message(paste("Simulation computing time:",round((t1-t0)/60,1),"min"))
message(paste("Total elapsed time:       ",round((t2-t0)/60,1),"min"))

