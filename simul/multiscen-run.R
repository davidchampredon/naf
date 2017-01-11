###
###   RUN MULTIPLE SCENARIOS AS DEFINED IN A FILE
###


# Load all default parameters:
source('utils-load-param.R')
source('utils-run.R')
source('utils-misc.R')

dir.results    <- dir.def('dir-def.csv')[['results']]
dir.save.rdata <- dir.def('dir-def.csv')[['rdata']]

# File defining parameters for various scenarios:
scen.list.file <- 'scenario-prm-list.csv'
x <- read.csv(scen.list.file)
x <- na.omit(x)
scen.id <- x$scenario_id

# Run the simulations on all scenarios:
for(i in seq_along(scen.id)){
	# Overwrite parameter values
	# associated with current scenario:
	ov <- overwrite.selected.param(filename = scen.list.file,
								   scen.id  = scenidx,
								   prm = prm, interv.prm = interv.prm)
	prm <- ov[['prm']]
	interv.prm <- ov[['interv.prm']]
	
	
	# Run the simulation for that scenario:
	run.simul(scen.id = scen.id[i], 
			  dir.save.rdata = dir.save.rdata,
			  force.light.output = TRUE)
}
