##################################################################
######
######    MINIMAL TEST FOR 'stiagent' LIBRARY
######
######
##################################################################

library(ggplot2)
library(plyr)
library(snowfall)

library(naf,lib.loc = "./lib")

source('analysis_tools.R')
source("../../SEmInR/SEmInR_deterministic.R")

sfInit(parallel = TRUE, cpu = 4)
sfLibrary(naf,lib.loc = "./lib")

t0 <- Sys.time()

prm <- list()
simul.prm <- list()

cr <- 1.0
doi_mean <- 3.0
R0 <- cr * doi_mean

prm[['debug_mode']] <- FALSE
prm[['dol_mean']] <- 2.0
prm[['doi_mean']] <- doi_mean

prm[['proba_move']] <- 0.999
prm[['homogeneous_contact']] <- TRUE
prm[['contact_rate']] <- cr

simul.prm[['horizon']] <- 70
simul.prm[['n_indiv']] <- 3E2
simul.prm[['initial_latent']] <- 2
simul.prm[['nt']] <- 4


# Deterministic SEIR
prm1 <- c(
	infectious_mean = prm[['doi_mean']],
	latent_mean     = prm[['dol_mean']],
	popSize         = simul.prm[['n_indiv']],
	R0              = R0)
prm2 <-  c(
	horizon     = simul.prm[['horizon']],
	nE = 1,
	nI = 1,
	init_I1 = simul.prm[['initial_latent']] ,
	n.time.steps = 500,
	per.capita = FALSE
)

### ODEs deterministic simulations:
sim.det <- simul.SEmInR.det(prm1, prm2)
df  <- sim.det$ts
det.prev <- df$Iall
det.t <- df$time

n.MC <- 10


loop.MC.naf <- function(iMC, prm, simul.prm) {
	prm[['rnd_seed']] <- 3*iMC
	return(naf_test_SEIR_vs_ODE(prm, simul.prm))
}
idx.apply <- 1:n.MC
res <- sfSapply(idx.apply, loop.MC.naf, prm, simul.prm,simplify = FALSE)
sfStop()

### merge MC iterations in one dataframe:
for(k in 1:n.MC){
	if(k==1) {
		sim.naf <- data.frame(res[[1]][['time_series']])
		sim.naf$mc <- 1
	}
	if(k>1) {
		tmp <- data.frame(res[[k]][['time_series']])
		tmp$mc <- k
		sim.naf <- rbind(sim.naf,tmp)
	}
}

CI <- 0.95
sim.naf.2 <- ddply(sim.naf,"time",summarize, 
				   prev = mean(prevalence),
				   prev.lo = quantile(prevalence,probs = 0.5-CI/2),
				   prev.hi = quantile(prevalence,probs = 0.5+CI/2))

### PLOTS ###

plot.population(data.frame(res[[1]][['population_final']]))
plot.epi.timeseries(data.frame(res[[1]][['time_series']]))

# compare naf vs ODEs
plot(det.t, det.prev, 
	 typ='l', col='red', 
	 lty = 2, lwd =2,
	 main = paste0("naf v.s. SEIR ODEs (nMC=",n.MC,", pop=",simul.prm[['n_indiv']] ,")"),
	 ylim=range(det.prev,sim.naf.2$prev.hi))
lines(sim.naf.2$time, sim.naf.2$prev)
lines(sim.naf.2$time, sim.naf.2$prev.lo, col='grey')
lines(sim.naf.2$time, sim.naf.2$prev.hi, col='grey')

legend(x = 'topleft',
	   legend = c('Deterministic','naf'),col = c(2,1),
	   lty=c(2,1))

t1 <- Sys.time()
message(paste("time elapsed:",round(t1-t0,1),"sec"))

