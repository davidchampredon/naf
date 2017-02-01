###
###   ANALYZE & COMPARE MULTIPLE SCENARIOS PREVIOUSLY RUN
###

library(tidyr)
library(plyr)
library(dplyr)
library(gridExtra)

source('utils-analyze.R')
source('utils-compare.R')
source('utils-misc.R')

dir.results    <- dir.def('dir-def.csv')[['results']]
dir.save.rdata <- dir.def('dir-def.csv')[['rdata']]


# Define 'scen.id' by retrieving 
# all successful simulations:
tmp     <- system(command = paste0('ls ',dir.save.rdata,'mc-simul-*.RData'),
				  intern = TRUE)
tmp     <- gsub(x = tmp, pattern = paste0(dir.save.rdata,'mc-simul-'),replacement = '')
scen.id.vec <- gsub(x = tmp, pattern = '.RData',   replacement = '')

bfirst <- TRUE

for(i in seq_along(scen.id.vec)){
	# In-depth analysis of the results:
	if(bfirst) analyze.simul.scen(scen.id = scen.id.vec[i],
								  dir.save.rdata = dir.save.rdata,
								  dir.results = dir.results)
	bfirst <- FALSE
	# Compare intervention of this scenario
	# with common baseline:
	compare.simul.scen(scen.id = scen.id.vec[i], 
					   dir.save.rdata = dir.save.rdata,
					   dir.results = dir.results)
}

# Merge results from all scenarios
# in file 'result-scen-all.RData':
res.all <- merge.result.scen(scen.id.vec, dir.save.rdata)

msg.ac <- 'Analysis and/or comparison of multi-scenarios done.'
print(msg.ac)
message(msg.ac)
