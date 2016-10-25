##################################################################
######
######    FIT AGE DISTRIBUTIONS GIVEN HOUSEHOLD SIZE
######
######    Note: This will DIRECTLY change distribution files
######    in ../data/households/
######
##################################################################

t0 <- as.numeric(Sys.time())

R.library.dir    <- '../Rlibrary/lib'
param.model.dir  <- '../param-model/'
data.dir         <- '../data/'

library(ggplot2)
library(plyr)
library(naf,lib.loc = R.library.dir)
library(snowfall)
library(nloptr)

source('../simul/analysis_tools.R')
source('../simul/build-world-det-sizegen.R')
source('../data/households/gen-ad-hhsz-fcts.R')
source('fit_fcts.R')
source(paste0(param.model.dir,'read-prm.R'))
source(paste0(param.model.dir,'read-world-prm.R'))
source(paste0(param.model.dir,'read-interv.R'))
source(paste0(param.model.dir,'read-sched.R'))

# Parameter file names
fname.prm.epi      <- 'prm-epi.csv'
fname.prm.simul    <- 'prm-simul.csv'
fname.prm.au       <- 'prm-au-ontario.csv'
fname.prm.interv.0 <- 'prm-interv-0.csv'
fname.prm.interv   <- 'prm-interv.csv'
fname.schedules    <- 'schedules.csv'
fname.prm.fit      <- 'prm-fit.csv'

prm        <- list()
simul.prm  <- list()
interv.prm <- list()
world.prm  <- list()
sched.prm  <- list()
fit.prm    <- list()

### ==== Model Parameters ====

prm[['debug_mode']] <- FALSE

prm       <- read.prm(filename = paste0(param.model.dir,fname.prm.epi), 
					  prm = prm)
simul.prm <- read.prm(filename = paste0(param.model.dir,fname.prm.simul), 
					  prm = simul.prm)

fit.prm   <- read.prm(filename = paste0(param.model.dir,fname.prm.fit), 
					  prm = fit.prm)


stoch_build_world <- simul.prm[['build_world_stoch']]

world.prm <- load.world.prm(filename = paste0(param.model.dir,fname.prm.au), 
							path.prm = param.model.dir)

world.prm[['id_region']]  <- 0
world.prm[['regionName']] <- "Canada"
world.prm[['unemployed_prop']] <- 0.10

sf           <- as.numeric(simul.prm[['scale_factor']])
world.prm    <- scale.world(1/sf, world.prm)
print(paste('World size reduced by',sf))

# age distributions, conditional on 
# household composition:
base.fname  <- paste0(data.dir,'households/hh_size_ad')
world.prm   <- c(world.prm,
				 load.age.hh(base.fname, 
				 			max.hh.size = max(world.prm[['max_hh_size']])) )
world.prm <- gen_world_ontario(world.prm)

# schedule time slices:
sched.prm[['timeslice']] <- c(1.0/24, 4.0/24, 4.0/24, 
							  1.0/24, 2.0/24, 12.0/24)
sched.prm[['sched_desc']]  <- read.schedules(paste0(param.model.dir,fname.schedules), return.names = FALSE)
sched.prm[['sched_names']] <- read.schedules(paste0(param.model.dir,fname.schedules), return.names = TRUE)

# Baseline, no vax intervention:
interv.prm <- load.interv.prm(paste0(param.model.dir,fname.prm.interv.0))


### ==== Constraint optimization ====

agemean.vec <- c(49, 51, 23.2, 59.4, 52.1, 20.2)     
a.vec       <- c(1.57, 7.07, 4.54, 3.27, 4.37, 0.95)
x0          <- c(agemean.vec, a.vec)
# WARNING: agemean boundaries MUST be 
# consistent with the hard-coded ones in
# function 'gen.all.ad.hhsz'
agemean.hi  <- c(69, 69, 69, 66,66,45)
agemean.lo  <- c(30, 19, 2, 19, 2, 2)
a.hi        <- 11
a.lo        <- 0.7
ub <- c( agemean.hi , rep(a.hi,length(agemean.hi)) )
lb <- c( agemean.lo , rep(a.lo,length(agemean.lo)) )

# Choose the optimizing method here:
optim.method <- 'NLOPT'       # 'ABC' , 'NLOPT'

if(optim.method=='ABC'){
abmfit <- optim_abc(n.abc = fit.prm[['fit_hhsz_age_ABC_n']], 
				   n.cpu  = fit.prm[['fit_hhsz_age_cpu']],
				   lb, 
				   ub,  
				   prm, 
				   simul.prm, 
				   interv.prm, 
				   world.prm, 
				   sched.prm)
xbest <- abmfit$x.best
}

if(optim.method=='NLOPT'){
xbest <- optim_nlopt(algo.name = 'NLOPT_GN_ISRES', 
					 x0, lb, ub,  
					 prm, 
					 simul.prm, 
					 interv.prm, 
					 world.prm, 
					 sched.prm,
					 fit.prm)
}
save.image(file='fit_hhsz_age.RData')

# Retrieve the optimal value
# and plot the 'best' distribution:
plot.best.fit(xbest,
			  prm        = prm, 
			  simul.prm  = simul.prm, 
			  interv.prm = interv.prm, 
			  world.prm  = world.prm, 
			  sched.prm  = sched.prm)

t1 <- as.numeric(Sys.time())
message(paste("Computing time:",round((t1-t0)/60,1),"min"))

