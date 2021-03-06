####
####   HELPER FUNCTIONS FOR FIT AGE DISTRIBUTION CONDIITINAL ON HOUSEHOLDS SIZE
####

# Polynomial regression
poly.reg <- function(x,y,deg) {
	df <- data.frame(x,y)
	L <- lm( y~ poly(x, deg, raw=TRUE), df)
	return(L$fitted.values)
}


# Error distance b/w 2 curves
err.dist <- function(x,y,p){
	n <- min(length(x),length(y))
	x <- unname(x[1:n])
	y <- unname(y[1:n])
	s <- 0
	for(i in 1:n) s <- s + (x[i]-y[i])^p
	return(s^(1/p))
}


# Snowfall wrapper:
run.snow.wrap <- function(seedMC,
						  prm, 
						  simul.prm, 
						  interv.prm, 
						  world.prm, 
						  sched.prm,
						  stoch_build_world){
	if(stoch_build_world){
		res <- naf_run(prm, 
					   simul.prm, 
					   interv.prm, 
					   world.prm, 
					   sched.prm,
					   seedMC)
	}
	else{
		res <- naf_run_det(prm, 
						   simul.prm, 
						   interv.prm, 
						   world.prm, 
						   sched.prm,
						   seedMC)
	}	
	return(res)
}

# Construct the simulation world and 
# return the simulated _global_ age distribution:
get.age.update.world <- function(prm, 
								 simul.prm, 
								 interv.prm, 
								 world.prm, 
								 sched.prm) 
{
	simul.prm[['build_world_only']] <- TRUE
	
	# Update world parameter
	# Relevant when age distributions
	# given household size previously updated:
	
	# Remove existing:
	idx <- which(grepl('pr_age_hh_',names(world.prm)))
	world.prm[idx] <- NULL
	# Update with new ones:
	world.prm   <- c(world.prm,
					 load.age.hh(base.fname, 
					 			max(world.prm[['max_hh_size']])) )
	world.prm <- gen_world_ontario(world.prm)
	
	# Note: put in a list format to be consistent
	# with what snowfall does. 
	# Here snowfall is not really usefull as we most likely
	# run with 1 MC (when deterministic built for world).
	# I had to delete snowfall in this function as ther was 
	# issues with embedded export of functions...
	res <- list(run.snow.wrap(1,
							  prm, 
							  simul.prm, 
							  interv.prm, 
							  world.prm, 
							  sched.prm,
							  stoch_build_world))
	
	# Average the global age distributions
	# across all MC iterations:
	a <- list()
	for(i in seq_along(res)){
		a[[i]] <- res[[i]][['ages']]
	}
	amx <- max(unlist(lapply(a, max)))+1
	sim.age <- 0:amx
	
	M <- matrix(nrow = length(sim.age)-1, ncol=length(res))
	for(i in seq_along(res)){
		h <- hist(a[[i]], breaks = sim.age, plot = F)
		sim.age.prop.i <- h$counts / sum(h$counts)	
		M[,i] <- sim.age.prop.i
	}
	sim.age.prop <- apply(M,1,FUN = mean)
	
	return(list(sim.age = sim.age, 
				sim.age.prop = sim.age.prop))
}


error.function <- function(agemean.vec,
						   a.vec,
						   prm, 
						   simul.prm, 
						   interv.prm, 
						   world.prm, 
						   sched.prm) {
	
	# Generate age (target) distributions given household size
	# with parameters 'agemean' and 'a':
	gen.all.ad.hhsz(agemean.vec, a.vec, 
					path.save='../data/households/')
	
	
	A <- get.age.update.world(prm        = prm, 
							  simul.prm  = simul.prm, 
							  interv.prm = interv.prm, 
							  world.prm  = world.prm, 
							  sched.prm  = sched.prm) 
	# Global age distribution:
	sim.age      <- A[['sim.age']]
	sim.age.prop <- A[['sim.age.prop']]
	sim.age2     <- sim.age[1:length(sim.age.prop)]
	
	# Retrieve the target global age distribution from data:
	target.ad <- read.csv('../data/ages/size-distrib-ages.csv')
	idx       <- which(target.ad$age <= max(sim.age))
	t.ad      <- target.ad[idx,]
	t.ad$prop <- t.ad$prop / sum(t.ad$prop)
	
	# Polynomial regression
	# (fit is made on this smooth regression
	#  rather than noisy simulation outputs & data)
	do.reg <- FALSE
	if(do.reg){
		deg      <- 7
		t.ad.reg <- poly.reg(x=t.ad$age, y=t.ad$prop, deg)
		sim.reg  <- poly.reg(x=sim.age2, y=sim.age.prop, deg)
	}
	if(!do.reg){
		t.ad.reg <- t.ad$prop
		sim.reg  <- sim.age.prop
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


eval_f <- function(x,
				   prm, 
				   simul.prm, 
				   interv.prm, 
				   world.prm, 
				   sched.prm){
	agemean.vec <- x[1:6]
	a.vec       <- x[7:12]
	EF <- error.function(agemean.vec, 
						 a.vec,
						 prm, 
						 simul.prm, 
						 interv.prm, 
						 world.prm, 
						 sched.prm)
	return(EF[['err']])
}


# Optimization using 'nlopt' library

optim_nlopt <- function(algo.name, x0, lb, ub,  
						prm, 
						simul.prm, 
						interv.prm, 
						world.prm, 
						sched.prm,
						fit.prm) {
	
	# for all options, see: nloptr.print.options()
	opts <- list('algorithm'   =algo.name, #'NLOPT_GN_ISRES', # NLOPT_GN_DIRECT_L, NLOPT_GN_ISRES,     NLOPT_GN_CRS2_LM', # 'NLOPT_LN_COBYLA',
				 'print_level' = 3,
				 'maxtime'     = 60*fit.prm[['fit_maxtime_minutes']])
	
	opt.val <- nloptr(x0 = x0, 
					  eval_f = eval_f, 
					  lb = lb, 
					  ub = ub, 
					  opts = opts,
					  prm        = prm, 
					  simul.prm  = simul.prm, 
					  interv.prm = interv.prm, 
					  world.prm  = world.prm, 
					  sched.prm  = sched.prm)
	return(opt.val$solution)
}

# Wrap for parallel execution:
snow.abc <- function(i,
					 MABC,
					 prm, 
					 simul.prm, 
					 interv.prm, 
					 world.prm, 
					 sched.prm){
	x <- MABC[i,]
	tryeval <- try(res <- eval_f(x = x,
								 prm        = prm, 
								 simul.prm  = simul.prm, 
								 interv.prm = interv.prm, 
								 world.prm  = world.prm, 
								 sched.prm  = sched.prm)
	)
	if(class(tryeval)=='try-error') res <- NA
	return(res)
}

optim_abc <- function(n.abc, 
					  n.cpu,
					  lb, ub,  
					  prm, 
					  simul.prm, 
					  interv.prm, 
					  world.prm, 
					  sched.prm){
	
	# Matrix of all parameters combinations:
	b <- length(lb) ; stopifnot(b==length(ub))
	MABC <- matrix(nrow=n.abc, ncol=b)
	for(j in 1:b){
		tmp <- seq(lb[j],ub[j],length.out = n.abc)
		MABC[,j] <- sample(tmp)
	}
	
	sfInit(parallel = (n.cpu>1), cpu = n.cpu)
	sfExportAll()
	sfLibrary(naf, lib.loc = R.library.dir) 
	res <- sfSapply(x          = 1:n.abc, 
					fun        = snow.abc,
					MABC       = MABC,
					prm        = prm, 
					simul.prm  = simul.prm, 
					interv.prm = interv.prm, 
					world.prm  = world.prm, 
					sched.prm  = sched.prm,
					simplify   = FALSE)
	sfStop()
	idx.best <- which.min(res)
	
	return(list(x.best = MABC[idx.best,], 
				fx.best = res[idx.best],
				x.all = MABC,
				fx.all = unlist(res)))
}


plot.best.fit <- function(xbest,
						  prm, 
						  simul.prm, 
						  interv.prm, 
						  world.prm, 
						  sched.prm) {
	
	EF <- error.function(agemean.vec = xbest[1:6],
						 a.vec       = xbest[7:12],
						 prm        = prm, 
						 simul.prm  = simul.prm, 
						 interv.prm = interv.prm, 
						 world.prm  = world.prm, 
						 sched.prm  = sched.prm)
	
	err          <- EF[['err']]
	t.ad         <- EF[['t.ad']]
	sim.age      <- EF[['sim.age']]
	sim.age.prop <- EF[['sim.age.prop']]
	target.reg   <- EF[['target.reg']]
	sim.reg      <- EF[['sim.reg']]
	
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
	
	try(read.all.ad.hhsz(path = '../data/households/'), silent=FALSE)
	dev.off()
}



plot.best.fit.manual <- function(x1,
								 x2,
								 prm, 
								 simul.prm, 
								 interv.prm, 
								 world.prm, 
								 sched.prm) {
	
	# First parameter set:
	
	EF <- error.function(agemean.vec = x1[1:6],
						 a.vec       = x1[7:12],
						 prm        = prm, 
						 simul.prm  = simul.prm, 
						 interv.prm = interv.prm, 
						 world.prm  = world.prm, 
						 sched.prm  = sched.prm)
	
	err          <- EF[['err']]
	t.ad         <- EF[['t.ad']]
	sim.age      <- EF[['sim.age']]
	sim.age.prop <- EF[['sim.age.prop']]
	target.reg   <- EF[['target.reg']]
	sim.reg      <- EF[['sim.reg']]
	
	
	# second parameter set:
	EF <- error.function(agemean.vec = x2[1:6],
						 a.vec       = x2[7:12],
						 prm        = prm, 
						 simul.prm  = simul.prm, 
						 interv.prm = interv.prm, 
						 world.prm  = world.prm, 
						 sched.prm  = sched.prm)
	
	err2          <- EF[['err']]
	t.ad2         <- EF[['t.ad']]
	sim.age2      <- EF[['sim.age']]
	sim.age.prop2 <- EF[['sim.age.prop']]
	target.reg2   <- EF[['target.reg']]
	sim.reg2      <- EF[['sim.reg']]
	
	
	par(mfrow=c(1,1),mar=rep(2,4))
	
	plot(t.ad$age, 
		 t.ad$prop, 
		 ylim= range(t.ad$prop,sim.age.prop,0),
		 main = paste('err =',round(err,5),'  err2 =',round(err2,5)),
		 cex=1.0, typ='p', col='lightgrey')
	lines(t.ad$age, target.reg, lty =2, lwd =2)
	grid()
	
	col1 <- 'pink'
	col2 <- 'brown'
	
	lines(sim.age,
		  sim.age.prop, 
		  col = col1,
		  pch=16, cex=0.7,  typ='o')
	lines(sim.age, sim.reg, lty =1, lwd=2, col=col1)
	lines(sim.age2,
		  sim.age.prop2, 
		  col = col2,
		  pch=16, cex=0.7,  typ='o')
	lines(sim.age2, sim.reg2, 
		  lty =1, lwd=2, col=col2)
	
	#try(read.all.ad.hhsz(path = '../data/households/'), silent=FALSE)
}




load.all.parameters <- function(R.library.dir,
								param.model.dir,
								data.dir,
								simul.dir) {
	
	library(naf,lib.loc = R.library.dir)
	
	source(paste0(simul.dir,'analysis_tools.R'))
	source(paste0(simul.dir,'build-world-det-sizegen.R'))
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
	
	prm        <- list()
	simul.prm  <- list()
	interv.prm <- list()
	world.prm  <- list()
	sched.prm  <- list()
	
	###  Model Parameters
	
	prm[['debug_mode']] <- FALSE
	
	prm       <- read.prm(filename = paste0(param.model.dir, fname.prm.epi), 
						  prm = prm)
	simul.prm <- read.prm(filename = paste0(param.model.dir, fname.prm.simul), 
						  prm = simul.prm)
	
	### World parameters
	
	world.prm <- load.world.prm(filename = paste0(param.model.dir, fname.prm.au), 
								path.prm = param.model.dir)
	world.prm[['id_region']]  <- 0
	world.prm[['regionName']] <- "Canada"
	
	# rescale world size:
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
	
	#  Schedules
	
	sched.prm[['timeslice']]   <- c(1.0/24, 4.0/24, 4.0/24, 
									1.0/24, 2.0/24, 12.0/24)
	sched.prm[['sched_desc']]  <- read.schedules(paste0(param.model.dir,fname.schedules), return.names = FALSE)
	sched.prm[['sched_names']] <- read.schedules(paste0(param.model.dir,fname.schedules), return.names = TRUE)
	
	### Intervention parameters
	
	# Baseline, no vax intervention:
	interv.prm.0 <- load.interv.prm(paste0(param.model.dir,fname.prm.interv.0))
	# Vax intervention:
	interv.prm   <- load.interv.prm(paste0(param.model.dir,fname.prm.interv))
	
	return(list(prm = prm,
				simul.prm = simul.prm,
				world.prm = world.prm,
				sched.prm = sched.prm,
				interv.prm.0 = interv.prm.0,
				interv.prm = interv.prm))
}




