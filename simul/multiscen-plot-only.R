### 
###   PLOTS SCENARIOS RESULTS
###

library(tidyr)
library(gridExtra)

source('utils-compare.R')
source('utils-misc.R')

dir.results    <- dir.def('dir-def.csv')[['results']]
dir.save.rdata <- dir.def('dir-def.csv')[['rdata']]

# ==== Load saved data ====

do.secondary <- FALSE

print('Loading saved RData...')
load(paste0(dir.save.rdata,'result-scen-all.RData'))
res.all <- list(main      = result.scen.all)

if(do.secondary) {
	load(paste0(dir.save.rdata,'secondary-all.RData'))
	res.all <- list(main      = result.scen.all, 
					secondary = secondary.all)
}
print('Saved RData loaded.')

# ==== Plots ====

print('Plotting main comparison...')


plot.rate.reduc(result.scen.all    = res.all[['main']], 
                dir                = dir.results,
                file.scen.prm.list = 'scenario-prm-list.csv')


figures.maintext(result.scen.all    = res.all[['main']], 
                 dir                = dir.results,
                 file.scen.prm.list = 'scenario-prm-list.csv')

if(FALSE){
    plot.multi.scen.res(result.scen.all    = res.all[['main']], 
					dir                = dir.results,
					file.scen.prm.list = 'scenario-prm-list.csv')
}

if(do.secondary){
	print('Plotting secondary comparison...')
	
	plot.secondary.res(res.all[['secondary']], 
					   dir = dir.results)
}



print('\n--> Scenario comparison plot completed.')
message('\n--> Scenario comparison plot completed.')
