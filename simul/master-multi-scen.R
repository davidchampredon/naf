##################################################################
######
######    MASTER SCRIPT TO RUN MULTIPLE SCENARIOS
######
######
##################################################################

# Load all default parameters:
source('utils-load-param.R')
source('utils-run.R')
source('utils-analyze.R')
source('utils-compare.R')

# File defining parameters for various scenarios:
scen.list.file <- 'scenario-prm-list.csv'
x <- read.csv(scen.list.file)
x <- na.omit(x)
scen.id <- x$scenario_id

# Run the simulations on all scenarios:
for(i in seq_along(scen.id)){
	# Overwrite parameter values
	# associated with current scenario:
	overwrite.selected.param(filename = scen.list.file, scen.id = scen.id[i]) 
	
	# Run the simulation for that scenario:
	run.simul(scen.id = scen.id[i])
	
	# In-depth analysis of the results:
	if(FALSE) analyze.simul.scen(scen.id = scen.id[i])
	
	# Compare intervention of this scenario
	# with common baseline:
	compare.simul.scen(scen.id = scen.id[i])
}

