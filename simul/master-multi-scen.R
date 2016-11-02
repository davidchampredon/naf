##################################################################
######
######    MASTER SCRIPT TO RUN MULTIPLE SCENARIOS
######
######
##################################################################

master0 <- as.numeric(Sys.time())

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
	overwrite.selected.param(filename = scen.list.file, 
							 scen.id = scen.id[i]) 
	
	# Run the simulation for that scenario:
	run.simul(scen.id = scen.id[i])
}

for(i in seq_along(scen.id)){

	# In-depth analysis of the results:
	if(FALSE) analyze.simul.scen(scen.id = scen.id[i])
	
	# Compare intervention of this scenario
	# with common baseline:
	compare.simul.scen(scen.id = scen.id[i])
}

# Merge results from all scenarios:
res.all <- merge.result.scen(scen.id)

# Plot results for all scenarios:
plot.multi.scen.res(res.all[['main']])
plot.secondary.res(res.all[['secondary']])

master1 <- as.numeric(Sys.time()) ; msgt <- paste("Time elapsed for master:",round((master1-master0)/60,1),'minutes.')
print(msgt)
message(msgt)
