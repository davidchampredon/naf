##################################################################
######
######    MINIMAL TEST FOR 'stiagent' LIBRARY
######
######
##################################################################

library(ggplot2)

library(naf,lib.loc = "./lib")

source('analysis_tools.R')
source("../../SEmInR/SEmInR_deterministic.R")




t0 <- Sys.time()

prm <- list()
simul.prm <- list()

cr <- 1.0
doi_mean <- 3.0
R0 <- cr * doi_mean

prm[['debug_mode']] <- TRUE
prm[['rnd_seed']] <- 123
prm[['dol_mean']] <- 2.0
prm[['doi_mean']] <- doi_mean

prm[['proba_move']] <- 0.999
prm[['homogeneous_contact']] <- TRUE
prm[['contact_rate']] <- cr

simul.prm[['horizon']] <- 70
simul.prm[['n_indiv']] <- 4000
simul.prm[['initial_latent']] <- 2



res <- naf_test(prm, simul.prm)

pop.final <- data.frame(res[['population_final']])
ts <- data.frame(res[['time_series']])


### PLOTS ###

plot.population(pop.final)
plot.epi.timeseries(ts)


t1 <- Sys.time()
message(paste("time elapsed:",round(t1-t0,1),"sec"))

