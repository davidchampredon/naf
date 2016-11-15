### 
###   PLOTS SCENARIOS RESULTS
###

source('utils-compare.R')
source('utils-misc.R')

dir.results    <- dir.def('dir-def.csv')[['results']]
dir.save.rdata <- dir.def('dir-def.csv')[['rdata']]

load(paste0(dir.save.rdata,'result-scen-all.RData'))
load(paste0(dir.save.rdata,'secondary-all.RData'))
res.all <- list(main      = result.scen.all, 
				secondary = secondary.all)

# Plot results for all scenarios:

plot.multi.scen.res(res.all[['main']], 
					dir = dir.results,
					file.scen.prm.list = 'scenario-prm-list.csv')

plot.secondary.res(res.all[['secondary']], 
				   dir = dir.results)

plot.rate.reduc(res.all[['main']], 
	  dir = dir.results,
	  file.scen.prm.list = 'scenario-prm-list.csv')


