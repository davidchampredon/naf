##################################################################
######
######    MINIMAL TEST FOR 'NAF' LIBRARY
######
######    
######
######
##################################################################

t0 <- as.numeric(Sys.time())

library(ggplot2)
library(plyr)

library(naf,lib.loc = "./lib")

source('analysis_tools.R')

save.plot.to.file <- T

prm <- list()
simul.prm <- list()
interv.prm <- list()
world.prm <- list()
sched.prm <- list()



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
simul.prm[['initial_latent']] <- 2
simul.prm[['popexport']] <- 1

###  ==== World parameters ====

world.prm[['id_au']] <- c(0,1)
world.prm[['name_au']] <- c("AUone", "AUtwo")
world.prm[['id_region']] <- 0
world.prm[['regionName']] <- "Region1"

world.prm[['n_hh']]     <- c(8, 7)
world.prm[['n_wrk']]    <- c(30, 33)
world.prm[['n_pubt']]   <- c(40, 30)
world.prm[['n_school']] <- c(20, 20)
world.prm[['n_hosp']]   <- c(1,1)
world.prm[['n_other']]  <- c(100, 50)

age.adults   <- c(22,33,44)
age.children <- c(11)
age.all      <- c(age.children, age.adults) 

world.prm[['pr_age_hh_00_val']] <- age.adults
world.prm[['pr_age_hh_10_val']] <- age.adults
world.prm[['pr_age_hh_11_val']] <- age.all
world.prm[['pr_age_hh_20_val']] <- age.adults
world.prm[['pr_age_hh_21_val']] <- age.adults
world.prm[['pr_age_hh_22_val']] <- age.all

world.prm[['pr_age_hh_00_proba']] <- c(0.2, 0.5, 0.3)
world.prm[['pr_age_hh_10_proba']] <- c(0.2, 0.5, 0.3)
world.prm[['pr_age_hh_11_proba']] <- c(0.2,0.1, 0.6, 0.1)
world.prm[['pr_age_hh_20_proba']] <- c(0.2, 0.5, 0.3)
world.prm[['pr_age_hh_21_proba']] <- c(0.2, 0.5, 0.3)
world.prm[['pr_age_hh_22_proba']] <- c(0.8,0.1, 0.05, 0.05)

world.prm[['hh_size']] <- c(1,2,3)
world.prm[['hh_size_proba']] <- c(0.3, 0.5, 0.2)
world.prm[['wrk_size']] <- c(5,20,40,60)
world.prm[['wrk_size_proba']] <- c(0.55, 0.3, 0.1, 0.05)
world.prm[['pubt_size']] <- c(20,50,60)
world.prm[['pubt_size_proba']] <- c(0.3,0.4,0.3)
world.prm[['school_size']] <- c(6, 8, 9)
world.prm[['school_size_proba']] <- c(0.7,0.2,0.1)
world.prm[['hosp_size']] <- c(10000)
world.prm[['hosp_size_proba']] <- c(1)
world.prm[['other_size']] <- c(10,50,100)
world.prm[['other_size_proba']] <- c(0.4,0.4,0.2)

sched.prm[['timeslice']] <- c(1.0/24, 4.0/24, 4.0/24, 1.0/24, 2.0/24, 12.0/24)


###  ==== Intervention parameters ====

interv.prm[['interv_name']] <- 'interv_test'
interv.prm[['interv_type']] <- 'vaccination'  # treatment cure vaccination
interv.prm[['interv_target']] <- 'susceptible'  # symptomatic  susceptible
interv.prm[['interv_start']] <- 15
interv.prm[['interv_end']] <- 999
interv.prm[['interv_cvg_rate']] <- 0.05
interv.prm[['interv_cvg_max_prop']] <- 0.9999


### ==== Run Simulation ====


res <- naf_run(prm, 
			   simul.prm, 
			   interv.prm, 
			   world.prm, 
			   sched.prm)

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



