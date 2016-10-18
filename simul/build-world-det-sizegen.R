

gen_world_ontario <- function(world.prm){
	### 
	### The goal of this function is to
	### construct a _STABLE_ world.
	### (by stable, I mean where the SP and 
	### individuals sizes are consistent, which
	### is not guaranteed when stochastically constructed)
	###
	
	set.seed(1234567)
	
	pathprm <- '../param-model/'
	
	get_real_data <- function(sptype) {
		X <- read.csv(paste0(pathprm,'prm-size-',sptype,'-ontario.csv'))
		X_val <- X$size
		X_proba <- X$prop
		X_n <- world.prm[[paste0('n_',sptype)]]
		return(list(val=X_val, proba=X_proba, n=X_n))
	}
	
	# Social places that are NOT explictly changed
	# from the real data:
	hh    <- get_real_data('hh')
	wrk   <- get_real_data('wrk')
	pubt  <- get_real_data('pubt')
	other <- get_real_data('other')
	hosp  <- get_real_data('hosp')
	
	size_hh    <- sample(x = hh[['val']], size = hh[['n']], replace = TRUE, prob = hh[['proba']])
	size_wrk   <- sample(x = wrk[['val']], size = wrk[['n']], replace = TRUE, prob = wrk[['proba']])
	size_pubt  <- sample(x = pubt[['val']], size = pubt[['n']], replace = TRUE, prob = pubt[['proba']])
	size_other <- sample(x = other[['val']], size = other[['n']], replace = TRUE, prob = other[['proba']])
	size_hosp  <- hosp[['val']] # doesn't work? ->> sample(x = hosp[['val']], size = hosp[['n']], replace = TRUE, prob = hosp[['proba']])
	
	# Social places that are explictly changed
	# from the real data:
	school <- get_real_data('school')
	size_school <- sample(x = school[['val']], size = school[['n']], replace = TRUE, prob = school[['proba']])
	
	# EXPLICIT CHANGE TO NUMBER AND/OR SIZES:
	size_school <- c(size_school, seq(500,2000,by=500))
	
	
	# Return the list of parameters:
	world.prm[['det_size_hh']]     <- list(size_hh)
	world.prm[['det_size_wrk']]    <- list(size_wrk)
	world.prm[['det_size_pubt']]   <- list(size_pubt)
	world.prm[['det_size_other']]  <- list(size_other)
	world.prm[['det_size_hosp']]   <- list(size_hosp)
	world.prm[['det_size_school']] <- list(size_school)
	
	return(world.prm)
}