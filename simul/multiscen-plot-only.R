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

print('Loading saved RData...')
load(paste0(dir.save.rdata,'result-scen-all.RData'))
res.all <- list(main = result.scen.all,
                main.ageGroup = result.scen.all.ageGroup)
print('Saved RData loaded.')

# ==== Plots ====

print('Plotting main comparison...')


plot.rate.reduc(df  = res.all[['main']], 
                dir = dir.results,
                file.scen.prm.list = 'scenario-prm-list.csv',
                do.ageGroup = FALSE)

plot.rate.reduc(df  = res.all[['main.ageGroup']], 
                dir = dir.results,
                file.scen.prm.list = 'scenario-prm-list.csv',
                do.ageGroup = TRUE)

try(expr = {
figures.maintext(result.scen.all    = res.all[['main']], 
                 dir                = dir.results,
                 file.scen.prm.list = 'scenario-prm-list.csv')
})


print('\n--> Scenario comparison plot completed.')
message('\n--> Scenario comparison plot completed.')
