##################################################################
######
######    MINIMAL TEST FOR 'NAF' LIBRARY
######
######    --> Does NAF simulate the right migration patterns in 
######    a simple 2-spatial location model?
######
######
##################################################################

library(ggplot2)
library(plyr)

library(naf,lib.loc = "./lib")

source('analysis_tools.R')

t0 <- Sys.time()

prm <- list()
simul.prm <- list()

cr <- 0.0
doi_mean <- 3.0
R0 <- cr * doi_mean

prm[['debug_mode']] <- FALSE
prm[['dol_mean']] <- 2.0
prm[['doi_mean']] <- doi_mean

prm[['proba_move']] <- 0.5
prm[['homogeneous_contact']] <- TRUE
prm[['contact_rate']] <- cr
prm[['rnd_seed']] <- 1234

simul.prm[['horizon']] <- 30
simul.prm[['n_indiv']] <- 50
simul.prm[['initial_latent']] <- 2
simul.prm[['nt']] <- 2


res <- naf_test_2_sp(prm, simul.prm)

ts_census <- as.data.frame(res[['time_series_census']])

pdf('plot_TEST_naf_2_sp.pdf')
g <- ggplot(ts_census)
g <- g + geom_step(aes(x=time, y=pop_present), size=1) 
g <- g + geom_step(aes(x=time, y=nS), colour='green4') 
g <- g + geom_step(aes(x=time, y=nE), colour='orange') 
g <- g + geom_step(aes(x=time, y=nIs), colour='red')
g <- g + geom_step(aes(x=time, y=nR), colour='blue') 
g <- g + facet_wrap(~id_sp)
plot(g)


t1 <- Sys.time()
message(paste("time elapsed:",round(t1-t0,1),"sec"))

dev.off()