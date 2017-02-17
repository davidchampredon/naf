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

file.scen.prm.list <- 'scenario-prm-list.csv'
dir <- dir.results

try(expr = {
    df  <- res.all[['main.ageGroup']]
    figures.maintext(df, dir, file.scen.prm.list)
    figure.tmp(result.scen.all, file.scen.prm.list)
})


if(FALSE){
    df  <- res.all[['main']]
    plot.rate.reduc(df, dir, file.scen.prm.list, do.ageGroup = FALSE)
    
    df  <- res.all[['main.ageGroup']]
    plot.rate.reduc(df, dir, file.scen.prm.list, do.ageGroup = TRUE)
}


print('\n--> Scenario comparison plot completed.')
message('\n--> Scenario comparison plot completed.')
