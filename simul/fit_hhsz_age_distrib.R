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

source('analysis_tools.R')
source(paste0(param.model.dir,'read-prm.R'))
source(paste0(param.model.dir,'read-world-prm.R'))
source(paste0(param.model.dir,'read-interv.R'))
source(paste0(param.model.dir,'read-sched.R'))

source('../data/households/gen-ad-hhsz-fcts.R')

# Parameter file names
fname.prm.epi      <- 'prm-epi.csv'
fname.prm.simul    <- 'prm-simul.csv'
fname.prm.au       <- 'prm-au-ontario.csv'
fname.prm.interv.0 <- 'prm-interv-0.csv'
fname.prm.interv   <- 'prm-interv.csv'
fname.schedules    <- 'schedules.csv'

prm        <- list()
simul.prm  <- list()
interv.prm <- list()
world.prm  <- list()
sched.prm  <- list()


### ==== Model Parameters ====

prm[['debug_mode']] <- FALSE

prm       <- read.prm(filename = paste0(param.model.dir,fname.prm.epi), 
					  prm = prm)
simul.prm <- read.prm(filename = paste0(param.model.dir,fname.prm.simul), 
					  prm = simul.prm)

###  ==== World parameters ====

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

max(world.prm[['max_hh_size']])

# schedule time slices:
sched.prm[['timeslice']] <- c(1.0/24, 4.0/24, 4.0/24, 
							  1.0/24, 2.0/24, 12.0/24)

sched.prm[['sched_desc']]  <- read.schedules(paste0(param.model.dir,fname.schedules), return.names = FALSE)
sched.prm[['sched_names']] <- read.schedules(paste0(param.model.dir,fname.schedules), return.names = TRUE)

###  ==== Intervention parameters ====

# Baseline, no vax intervention:
interv.prm <- load.interv.prm(paste0(param.model.dir,fname.prm.interv.0))

# Snowfall wrapper:
run.snow.wrap <- function(seedMC,
						  prm, 
						  simul.prm, 
						  interv.prm, 
						  world.prm, 
						  sched.prm){
	res <- naf_run(prm, 
				   simul.prm, 
				   interv.prm, 
				   world.prm, 
				   sched.prm,
				   seedMC)	
	return(res)
}


### ==== Run Simulation ====

error.function <- function(agemean.vec,a.vec) {
	
	# Generate age (target) distributions given household size
	# with parameters 'agemean' and 'a':
	gen.all.ad.hhsz(agemean.vec, a.vec, 
					path='../data/households/')
	# Construct the simulation world and 
	# return the simulated _global_ age distribution:
	get.age.update.world <- function(prm, 
									 simul.prm, 
									 interv.prm, 
									 world.prm, 
									 sched.prm) {
		n.MC  <- simul.prm[['fit_hhsz_age_mc']]
		n.cpu <- simul.prm[['fit_hhsz_age_cpu']]
		seeds <- 1:n.MC
		
		# Update world parameter
		# Relevant when age distributions
		# given household size previously updated:
		world.prm   <- c(world.prm,
						 load.age.hh(base.fname, 
						 			max(world.prm[['max_hh_size']])) )
		
		sfInit(parallel = (n.cpu>1), cpu = n.cpu)
		sfLibrary(naf,lib.loc = R.library.dir) 
		
		simul.prm[['build_world_only']] <- TRUE
		
		res <- sfSapply(x          = seeds, 
						fun        = run.snow.wrap,
						prm        = prm, 
						simul.prm  = simul.prm, 
						interv.prm = interv.prm, 
						world.prm  = world.prm, 
						sched.prm  = sched.prm,
						simplify   = FALSE)
		
		sfStop()
		
		# Average the global age distributions
		# across all MC iterations:
		a <- list()
		for(i in seq_along(res)){
			a[[i]] <- res[[i]][['ages']]
		}
		amx <- max(unlist(lapply(a, max)))
		sim.age <- 1:amx
		M <- matrix(nrow = amx-1, ncol=length(res))
		for(i in seq_along(res)){
			h <- hist(a[[i]], breaks = sim.age, plot = F)
			sim.age.prop.i <- h$counts / sum(h$counts)	
			M[,i] <- sim.age.prop.i
		}
		sim.age.prop <- apply(M,1,FUN = mean)
		
		return(list(sim.age = sim.age, 
					sim.age.prop = sim.age.prop))
	}
	A <- get.age.update.world(prm, 
							  simul.prm, 
							  interv.prm, 
							  world.prm, 
							  sched.prm) 
	# Global age distribution:
	sim.age      <- A[['sim.age']]
	sim.age.prop <- A[['sim.age.prop']]
	sim.age2     <- sim.age[1:length(sim.age.prop)]
	
	# Polynomial regression
	poly.reg <- function(x,y,deg) {
		df <- data.frame(x,y)
		L <- lm( y~ poly(x, deg, raw=TRUE), df)
		return(L$fitted.values)
	}
	
	# Retrieve the target global age distribution from data:
	target.ad <- read.csv('../data/ages/size-distrib-ages.csv')
	idx       <- which(target.ad$age <= max(sim.age))
	t.ad      <- target.ad[idx,]
	t.ad$prop <- t.ad$prop / sum(t.ad$prop)
	
	# Polynomial regression
	# (fit is made on this smooth regression
	#  rather than noisy simulation outputs & data)
	deg      <- 7
	t.ad.reg <- poly.reg(x=t.ad$age, y=t.ad$prop, deg)
	sim.reg  <- poly.reg(x=sim.age2, y=sim.age.prop, deg)
	
	err.dist <- function(x,y,p){
		n <- min(length(x),length(y))
		x <- unname(x[1:n])
		y <- unname(y[1:n])
		s <- 0
		for(i in 1:n) s <- s + (x[i]-y[i])^p
		return(s^(1/p))
	}
	
	# Error value (to be minimized):
	err <- err.dist(sim.reg, t.ad.reg, p=2)
	
	return(list(err = err,
				t.ad = t.ad,
				sim.age = sim.age2,
				sim.age.prop = sim.age.prop,
				target.reg = t.ad.reg,
				sim.reg = sim.reg))
}

agemean.vec <- c(65,59,59,35,35,12)   #c(60,50,50,35,35,10)
a.vec       <- c(1.5,2.5,2.5,1.8,1.8,2.1)

EF <- error.function(agemean.vec,a.vec)

err <- EF[['err']]

t.ad <- EF[['t.ad']]
sim.age <- EF[['sim.age']]
sim.age.prop <- EF[['sim.age.prop']]
target.reg <- EF[['target.reg']]
sim.reg <- EF[['sim.reg']]


### ==== Constraint optimization ====

eval_f <- function(x){
	agemean.vec <- x[1:6]
	a.vec       <- x[7:12]
	EF <- error.function(agemean.vec,a.vec)
	return(EF[['err']])
}

agemean.vec <- c(65,59,59,35,35,12)   #c(60,50,50,35,35,10)
a.vec       <- c(1.5,2.5,2.5,1.8,1.8,2.1)

x0 <- c(agemean.vec, a.vec)
lb <- c( c(40,30,30,25,25,5) ,  rep(0.1,6) )
ub <- c( c(70,70,70,65,65,45) , rep(10,6) )

opts <- list('algorithm'   ='NLOPT_GN_CRS2_LM', # 'NLOPT_LN_COBYLA',
			 'print_level' = 1,
			 'maxtime'     = 60*simul.prm[['fit_maxtime_minutes']])
# nloptr.print.options()

opt.val <- nloptr(x0 = x0,eval_f = eval_f, lb = lb, ub=ub, opts = opts)

xx <- opt.val$solution


# Retrieve the optimal value
# and plot the 'best' distribution:
EF <- error.function(agemean.vec = xx[1:6],
					 a.vec = xx[7:12])

err <- EF[['err']]

t.ad <- EF[['t.ad']]
sim.age <- EF[['sim.age']]
sim.age.prop <- EF[['sim.age.prop']]
target.reg <- EF[['target.reg']]
sim.reg <- EF[['sim.reg']]

pdf('plot-fit-hhsz-age.pdf', width = 15, height = 10)
par(mfrow=c(1,1))
plot(t.ad$age, t.ad$prop, 
	 ylim= range(t.ad$prop,sim.age.prop,0),
	 main = paste('err =',round(err,5)),
	 cex=1.0, typ='p', col='lightgrey')
lines(t.ad$age, target.reg, lty =2, lwd =2)
grid()
lines(sim.age,
	  sim.age.prop, 
	  col = 'pink',
	  pch=16, cex=0.7,  typ='o')
lines(sim.age, sim.reg, lty =1, lwd=2, col='red')

read.all.ad.hhsz(path = '../data/households/')
dev.off()

save.image(file='fit_hhsz_age.RData')
t1 <- as.numeric(Sys.time())
message(paste("Computing time:",round((t1-t0)/60,1),"min"))


