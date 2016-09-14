##################################################################
######
######    TEST FOR 'NAF' LIBRARY
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

prm[['dol_distrib']] <- "lognorm"
prm[['doi_distrib']] <- "lognorm"
prm[['doh_distrib']] <- "exp"

prm[['dol_mean']] <- 2.0
prm[['doi_mean']] <- doi_mean
prm[['doh_mean']] <- 5.0

prm[['dol_var']] <- 2.1
prm[['doi_var']] <- 2.5
prm[['doh_var']] <- 1.00

prm[['proba_move']] <- 1
prm[['proba_change_sp_other']] <- 0.11

prm[['proba_death_prm_1']] <- 0.08
prm[['proba_death_prm_2']] <- 0.75
prm[['proba_death_prm_3']] <- 9.999

prm[['frailty_0']]        <- 0.60
prm[['frailty_min']]      <- 0.15
prm[['frailty_agemin']]   <- 30.0
prm[['frailty_agepivot']] <- 60
prm[['frailty_pivot']]    <- 0.50
prm[['frailty_sd']]       <- 0.1
prm[['frailty_powerChild']] <- 3


prm[['homogeneous_contact']] <- FALSE
prm[['contact_rate_mean']]   <- cr
prm[['contact_rate_stddev']] <- 0.4

prm[['contact_ratio_age_1_10']]     <- 2.0
prm[['contact_ratio_age_10_16']]    <- 1.5
prm[['contact_ratio_age_over_65']]  <- 0.8

prm[['contact_ratio_sp_household']]     <- 2.0
prm[['contact_ratio_sp_pubTransport']]  <- 1.75


prm[['asymptom_infectiousness_ratio']] <- 0.8
prm[['treat_doi_reduc']] <- 1.1123
prm[['treat_reduc_infect_mean']] <- 0.1

prm[['vax_imm_incr']]   <- 0.4
prm[['vax_frail_incr']] <- 0.2
prm[['vax_lag_full_efficacy']] <- 12


### ==== Simulation parameters ====

simul.prm[['rnd_seed']] <- 1234561
simul.prm[['horizon']] <- 300
simul.prm[['initial_latent']] <- 8
simul.prm[['popexport']] <- 1

###  ==== World parameters ====

world.prm[['id_au']] <- c(0,1)
world.prm[['name_au']] <- c("AUone", "AUtwo")
world.prm[['id_region']] <- 0
world.prm[['regionName']] <- "Region1"

mult <- 4

# Number of social places in each area unit:
world.prm[['n_hh']]     <- c(1200, 1100) * mult
world.prm[['n_wrk']]    <- c(300, 330) * mult
world.prm[['n_pubt']]   <- c(400, 300) * mult
world.prm[['n_school']] <- c(200, 200) * mult
world.prm[['n_hosp']]   <- c(1,1)
world.prm[['n_other']]  <- c(1000, 500) * mult

# define coarse age groups and their distributions:
age.adult    <- seq(19,98,by=1)
age.children <- seq(1,18,by=1)
age.all      <- c(age.children, age.adult) 

n.age.adult <- length(age.adult)
n.age.all   <- length(age.all)

p.adult     <- synthetic_age_adult(age.adult)   
p.all       <- c(rep(p.adult[1],18), p.adult) 
p.all       <- p.all/sum(p.all)
p.all.child <- rep(0, n.age.all)
p.all.child[1:20] <- 1/20
p.all.child <- p.all.child/sum(p.all.child)

if(FALSE){
	plot(age.all, p.all, typ='l',
		 lwd = 2,
		 ylim=c(0,max(p.all,p.adult,p.all.child)))
	lines(age.adult, p.adult, lwd = 3,lty = 2, col = 'red')
	lines(age.all, p.all.child, lwd = 3,lty = 3)
}

# age distributions, conditional on 
# household composition:
world.prm[['pr_age_hh_00_val']] <- age.adult   
world.prm[['pr_age_hh_10_val']] <- age.adult   
world.prm[['pr_age_hh_11_val']] <- age.all     
world.prm[['pr_age_hh_20_val']] <- age.adult
world.prm[['pr_age_hh_21_val']] <- age.adult
world.prm[['pr_age_hh_22_val']] <- age.all

world.prm[['pr_age_hh_00_proba']] <- p.adult    # <-- age distribution of household of size 1
world.prm[['pr_age_hh_10_proba']] <- p.adult    # <-- age distribution of oldest indiv in household of size 2
world.prm[['pr_age_hh_11_proba']] <- p.all      # <-- age distribution of 2nd oldest indiv in household of size 2
world.prm[['pr_age_hh_20_proba']] <- p.adult
world.prm[['pr_age_hh_21_proba']] <- p.adult
world.prm[['pr_age_hh_22_proba']] <- p.all.child

# size distribution of all social places:
world.prm[['hh_size']] <- c(1,2,3)
world.prm[['hh_size_proba']] <- c(0.2, 0.4, 0.4)
world.prm[['wrk_size']] <- c(5,20,40,60)
world.prm[['wrk_size_proba']] <- c(0.55, 0.3, 0.1, 0.05)
world.prm[['pubt_size']] <- c(30,60,120)
world.prm[['pubt_size_proba']] <- c(0.3,0.4,0.3)
world.prm[['school_size']] <- c(100, 200, 300)
world.prm[['school_size_proba']] <- c(0.7,0.2,0.1)
world.prm[['hosp_size']] <- c(50000)
world.prm[['hosp_size_proba']] <- c(1)
world.prm[['other_size']] <- c(100,200,300)
world.prm[['other_size_proba']] <- c(0.4,0.4,0.2)

# schedule time slices:
sched.prm[['timeslice']] <- c(1.0/24, 4.0/24, 4.0/24, 1.0/24, 2.0/24, 12.0/24)


###  ==== Intervention parameters ====

interv.prm[['interv_name']]         <- 'interv_test'
interv.prm[['interv_type']]         <- 'vaccination'  # treatment cure vaccination
interv.prm[['interv_target']]       <- 'susceptible'  # symptomatic  susceptible
interv.prm[['interv_start']]        <- 20
interv.prm[['interv_end']]          <- 999
interv.prm[['interv_cvg_rate']]     <- 0.05
interv.prm[['interv_cvg_max_prop']] <- 0.9999


### ==== Run Simulation ====

try(
	res <- naf_run(prm, 
				   simul.prm, 
				   interv.prm, 
				   world.prm, 
				   sched.prm)
)

t1 <- as.numeric(Sys.time())

# ==== Process results ====

message("Processing results...")

ts <- as.data.frame(res[['time_series']])

# JUST ONE SOCIAL PLACE--->   pop <- as.data.frame(res[['population_final']])
world0 <- res[['world']]
z <- lapply(world0, as.data.frame)
pop <- do.call('rbind',z)

ws <- ddply(pop, c('id_au','sp_type'), summarize, 
			n_sp = length(id_sp),
			n_indiv = length(id_indiv))


### ==== PLOTS ==== 

message("Plotting results...")

if (save.plot.to.file) pdf('plot_TEST_naf.pdf', width = 30,height = 20)

try( plot.population(pop),  silent = T)
try( plot.n.contacts(res$track_n_contacts),  silent = T)
try( plot.epi.timeseries(ts),  silent = T)
try( grid.arrange(plot.ts.sp(res$time_series_sp),
				 plot.ts.sp(res$time_series_sp, facets = T)), 
	 silent = T)
try( plot.age.contact.matrix(res$wiw_ages),  silent = T)

if (save.plot.to.file) dev.off()


# End
t2 <- as.numeric(Sys.time())
message(paste("Simulation computing time:",round((t1-t0)/60,1),"min"))
message(paste("Total elapsed time:       ",round((t2-t0)/60,1),"min"))

