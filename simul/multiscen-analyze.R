###
###   ANALYZE & COMPARE MULTIPLE SCENARIOS PREVIOUSLY RUN
###

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
scen.id <- gsub(x = tmp, pattern = '.RData',   replacement = '')

for(i in seq_along(scen.id)){
	# In-depth analysis of the results:
	if(FALSE) analyze.simul.scen(scen.id = scen.id[i],
								 dir.save.rdata = dir.save.rdata)
	# Compare intervention of this scenario
	# with common baseline:
	compare.simul.scen(scen.id = scen.id[i], 
					   dir.save.rdata = dir.save.rdata,
					   dir.results = dir.results)
}

# Merge results from all scenarios
# in file 'result-scen-all.RData':
res.all <- merge.result.scen(scen.id, dir.save.rdata)

# Plot results for all scenarios:
plot.multi.scen.res(res.all[['main']], 
					dir=dir.results,
					file.scen.prm.list = 'scenario-prm-list.csv')

plot.secondary.res(res.all[['secondary']], 
				   dir = dir.results)

plot2(res.all[['main']], 
	  dir = dir.results,
	  file.scen.prm.list = 'scenario-prm-list.csv')

msg.ac <- 'Analysis and/or comparison of multi-scenarios done.'
print(msg.ac)
message(msg.ac)
