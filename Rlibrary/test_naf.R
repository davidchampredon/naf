##################################################################
######
######    MINIMAL TEST FOR 'NAF' LIBRARY
######
######    --> Does NAF simulate correctly hospitalization?
######
######
##################################################################

library(ggplot2)
library(plyr)

library(naf,lib.loc = "./lib")

source('analysis_tools.R')

save.plot.to.file <- T

t0 <- as.numeric(Sys.time())



### ==== Model Parameters ====

cr <- 2.1
doi_mean <- 4
R0 <- cr * doi_mean

prm[['debug_mode']] <- F

prm[['dol_distrib']] <- "exp"
prm[['doi_distrib']] <- "exp"
prm[['doh_distrib']] <- "exp"

prm[['dol_mean']] <- 2.0
prm[['doi_mean']] <- doi_mean
prm[['doh_mean']] <- 5.0

prm[['proba_move']] <- 1

prm[['proba_death_prm_1']] <- 0.08
prm[['proba_death_prm_2']] <- 0.75
prm[['proba_death_prm_3']] <- 9.999


prm[['homogeneous_contact']] <- FALSE
prm[['contact_rate']] <- cr
prm[['asymptom_infectiousness_ratio']] <- 0.8
prm[['treat_doi_reduc']] <- 1.1123
prm[['treat_reduc_infect_mean']] <- 0.1

prm[['vax_imm_incr']]   <- 0.4
prm[['vax_frail_incr']] <- 0.2
prm[['vax_lag_full_efficacy']] <- 12


### ==== Simulation parameters ====

simul.prm[['rnd_seed']] <- 1234
simul.prm[['horizon']] <- 300
simul.prm[['n_indiv']] <- 5000
simul.prm[['initial_latent']] <- 2
simul.prm[['nt']] <- 2
simul.prm[['popexport']] <- 1

# intervention:

simul.prm[['interv_name']] <- 'interv_test'
simul.prm[['interv_type']] <- 'vaccination'  # treatment cure vaccination
simul.prm[['interv_target']] <- 'susceptible'  # symptomatic  susceptible
simul.prm[['interv_start']] <- 15
simul.prm[['interv_end']] <- 999
simul.prm[['interv_cvg_rate']] <- 0.05
simul.prm[['interv_cvg_max_prop']] <- 0.9999


### ==== Run Simulation ====

res <- naf_test(prm, simul.prm)

ts <- as.data.frame(res[['time_series']])
ts_census <- as.data.frame(res[['time_series_census']])
pop <- as.data.frame(res[['population_final']])


### ==== PLOTS ==== 

if (save.plot.to.file) pdf('plot_TEST_naf.pdf', width = 30,height = 20)

plot.population(pop)

plot.epi.timeseries(ts)

g <- ggplot(ts_census)
g <- g + geom_point(aes(x=time, y=pop_present), size=1) 
g <- g + geom_step(aes(x=time, y=nS), colour='green4') 
g <- g + geom_step(aes(x=time, y=nE), colour='orange') 
g <- g + geom_step(aes(x=time, y=nIs), colour='red')
g <- g + geom_step(aes(x=time, y=nR), colour='blue') 
g <- g + facet_wrap(~id_sp)
plot(g)




t1 <- as.numeric(Sys.time())
message(paste("time elapsed:",round(t1-t0,1),"sec"))

if (save.plot.to.file) dev.off()