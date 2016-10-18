##################################################################
######
######    BUILD THE WORLD DETERMINISTICALLY
######
######    
######
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

source('analysis_tools.R')
source('build-world-det-sizegen.R')
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

do.plot           <- F# TRUE 
save.plot.to.file <- TRUE

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
world.prm   <- c(world.prm,
				 load.age.hh(base.fname, 
				 			max.hh.size = max(world.prm[['max_hh_size']])) )

world.prm <- gen_world_ontario(world.prm)


# schedule time slices:

sched.prm[['timeslice']] <- c(1.0/24, 4.0/24, 4.0/24, 
							  1.0/24, 2.0/24, 12.0/24)

sched.prm[['sched_desc']]  <- read.schedules(paste0(param.model.dir,fname.schedules), return.names = FALSE)
sched.prm[['sched_names']] <- read.schedules(paste0(param.model.dir,fname.schedules), return.names = TRUE)

###  ==== Intervention parameters ====

# Baseline, no vax intervention:
interv.prm.0 <- load.interv.prm(paste0(param.model.dir,fname.prm.interv.0))
# Vax intervention:
interv.prm   <- load.interv.prm(paste0(param.model.dir,fname.prm.interv))


### === Snowfall wrapper ===

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


### ==== BUILD WORLD ====

stoch_build_world <- FALSE
simul.prm[['build_world_only']] <- 1  # Must be 1, no simulation
n.MC  <- 1# simul.prm[['mc']]
n.cpu <- 1 #simul.prm[['cpu']]
seeds <- 1:n.MC

sfInit(parallel = (n.cpu>1), cpu = n.cpu)
sfLibrary(naf,lib.loc = R.library.dir) 

x <- sfSapply(seeds, run.snow.wrap,
			  prm        = prm, 
			  simul.prm  = simul.prm, 
			  interv.prm = interv.prm.0, 
			  world.prm  = world.prm, 
			  sched.prm  = sched.prm,
			  stoch_build_world = stoch_build_world,
			  simplify   = FALSE)

sfStop()

save.image(file='buildworld.RData')


### ==== Display world ====

sp.to.string <- function(i) {
	# WARNING: order matters!
	# MUST be same order as enum SPType definition (C++).
	if (i==0) res = 'household';
	if (i==1) res = 'workplace'; 
	if (i==2) res = 'school';
	if (i==3) res = 'other';
	if (i==4) res = 'hospital';
	if (i==5) res = 'pubTransp';
	# // [add here newly defined SPs...]
	# // if(i==6) res = SP_xxx;
	return(res)
}

sp      <- as.data.frame(x[[1]]$census_sp)
sp$type <- sapply(sp$sp_type,FUN = sp.to.string) 
sp.cnt  <- ddply(sp,c('type'),summarize, cnt=length(sp_id))

sp.cnt$target <- NA
sp.cnt$target[sp.cnt$type=='household'] <- world.prm[['n_hh']]
sp.cnt$target[sp.cnt$type=='workplace'] <- world.prm[['n_wrk']]
sp.cnt$target[sp.cnt$type=='other']     <- world.prm[['n_other']]
sp.cnt$target[sp.cnt$type=='school']    <- world.prm[['n_school']]
sp.cnt$target[sp.cnt$type=='pubTransp'] <- world.prm[['n_pubt']]
sp.cnt$target[sp.cnt$type=='hospital']  <- world.prm[['n_hosp']]
sp.cnt$prop.filled <- sp.cnt$cnt/sp.cnt$target
sp.cnt

g <- ggplot(sp.cnt) + geom_bar(aes(x=type, y=prop.filled), stat='identity')
g <- g + ggtitle('Proportion filled w.r.t. target sizes')

sp.dsz  <- ddply(sp,c('type','nlinked'), summarize, n=length(nlinked))

q <- join(sp.dsz, sp.cnt,by='type')
q$prop <- q$n/q$cnt
q$source <- 'simulated'

wp <- data.frame(type = 'household', nlinked=world.prm[['hh_size']], prop=world.prm[['hh_size_proba']])
wp <- rbind(wp, 
			data.frame(type = 'workplace', nlinked=world.prm[['wrk_size']], prop=world.prm[['wrk_size_proba']]))
wp <- rbind(wp, 
			data.frame(type = 'school', nlinked=world.prm[['school_size']], prop=world.prm[['school_size_proba']]))
wp <- rbind(wp, 
			data.frame(type = 'pubTransp', nlinked=world.prm[['pubt_size']], prop=world.prm[['pubt_size_proba']]))
wp <- rbind(wp, 
			data.frame(type = 'other', nlinked=world.prm[['other_size']], prop=world.prm[['other_size_proba']]))

wp$type <- as.character(wp$type)
wp$source <- 'target'

df <- rbind(q[,names(wp)],wp)
df$lt <- 1
df$lt[df$source=='target'] <- 2
g2 <- ggplot(df,aes(x=nlinked, y=prop, colour=source, shape=source)) + geom_point() + geom_line(aes(linetype=factor(lt))) + facet_wrap(~type, scale='free')
g2 <- g2 + ggtitle('Size Distribution: simulated v.s. target')

pdf('plot-build.pdf', width = 10, height = 8)
plot(g)
plot(g2)
dev.off()


t1 <- as.numeric(Sys.time())
message(paste("Computing time:",round((t1-t0)/60,1),"min"))


