##################################################################
######
######    TEST FOR 'NAF' LIBRARY
######
######    
######
######
##################################################################

t0 <- as.numeric(Sys.time())

R.library.dir    <- '../Rlibrary/lib'
param.model.dir  <- '../param-model/'

library(ggplot2)
library(plyr)
library(naf,lib.loc = R.library.dir)

source('analysis_tools.R')
source(paste0(param.model.dir,'read-prm.R'))

save.plot.to.file <- T

prm <- list()
simul.prm <- list()
interv.prm <- list()
world.prm <- list()
sched.prm <- list()


### ==== Model Parameters ====

cr <- 2.1

prm[['debug_mode']] <- F
prm       <- read.prm(filename = paste0(param.model.dir,'prm-epi.csv'), prm = prm)
simul.prm <- read.prm(filename = paste0(param.model.dir,'prm-simul.csv'), prm = simul.prm)


###  ==== World parameters ====

world.prm[['id_au']] <- c(0,1)
world.prm[['name_au']] <- c("AUone", "AUtwo")
world.prm[['id_region']] <- 0
world.prm[['regionName']] <- "Region1"

mult <- 2

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

do.plot <- FALSE 

if(do.plot){
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
}

# End
t2 <- as.numeric(Sys.time())
message(paste("Simulation computing time:",round((t1-t0)/60,1),"min"))
message(paste("Total elapsed time:       ",round((t2-t0)/60,1),"min"))

