###
###   RUN ONE SCENARIOS AS DEFINED IN A FILE
###

# >>> DELETE WHEN SURE <<<

# run.one.scenario <- function(scenidx) {
# 	
# 	# Load all default parameters:
# 	source('utils-load-param.R')
# 	source('utils-run.R')
# 	source('utils-misc.R')
# 	
# 	dir.results    <- dir.def('dir-def.csv')[['results']]
# 	dir.save.rdata <- dir.def('dir-def.csv')[['rdata']]
# 	
# 	# File defining parameters for various scenarios:
# 	scen.list.file <- 'scenario-prm-list.csv'
# 
# 	# Overwrite parameter values
# 	# associated with current scenario:
# 	ov <- overwrite.selected.param(filename   = scen.list.file,
# 								   scen.id    = scenidx,
# 								   prm        = prm, 
# 								   interv.prm = interv.prm)
# 	prm        <- ov[['prm']]
# 	interv.prm <- ov[['interv.prm']]
# 
# 	# Run the simulation for that scenario:
# 	run.simul(scen.id = scenidx, 
# 			  dir.save.rdata = dir.save.rdata,
# 			  force.light.output = TRUE)
# }
# 
