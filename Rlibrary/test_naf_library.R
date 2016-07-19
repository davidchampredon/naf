##################################################################
######
######    MINIMAL TEST FOR 'stiagent' LIBRARY
######
######
##################################################################

library(ggplot2)

library(naf,lib.loc = "./lib")

source('analysis_tools.R')

t0 <- Sys.time()

prm <- list()
prm[['debug_mode']] <- TRUE
prm[['rnd_seed']] <- 123
prm[['dol_mean']] <- 2.0
prm[['doi_mean']] <- 3.0
prm[['horizon']] <- 50
prm[['n_indiv']] <- 1000
prm[['proba_move']] <- 0.999
prm[['homogeneous_contact']] <- TRUE
prm[['contact_rate']] <- 0.08

res <- naf_test(prm)

pop <- data.frame(res[['population_final']])
ts <- data.frame(res[['time_series']])


### PLOTS ###

plot.epi.timeseries(ts)


t1 <- Sys.time()
message(paste("time elapsed:",round(t1-t0,1),"sec"))

