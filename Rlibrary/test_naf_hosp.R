##################################################################
######
######    MINIMAL TEST FOR 'NAF' LIBRARY
######
######    --> Does NAF simulate correctly hospitalization?
######
######
##################################################################

library(ggplot2)
library(plyr)

library(naf,lib.loc = "./lib")

source('analysis_tools.R')

save.plot.to.file <- F

t0 <- as.numeric(Sys.time())

prm <- list()
simul.prm <- list()

cr <- 2.9
doi_mean <- 6.0
R0 <- cr * doi_mean

prm[['debug_mode']] <- F

prm[['dol_distrib']] <- "exp"
prm[['doi_distrib']] <- "exp"
prm[['doh_distrib']] <- "exp"

prm[['dol_mean']] <- 2.0
prm[['doi_mean']] <- doi_mean
prm[['doh_mean']] <- 5.0

prm[['proba_move']] <- 1
prm[['homogeneous_contact']] <- FALSE
prm[['contact_rate']] <- cr
prm[['rnd_seed']] <- 12345

simul.prm[['horizon']] <- 40
simul.prm[['n_indiv']] <- 1000
simul.prm[['initial_latent']] <- 2
simul.prm[['nt']] <- 2
simul.prm[['popexport']] <- 1

res <- naf_test_hosp(prm, simul.prm)

ts_census <- as.data.frame(res[['time_series_census']])


### PLOTS ###

if (save.plot.to.file) pdf('plot_TEST_naf_hosp.pdf')

popf <- as.data.frame(res[['population_final']])
plot.population(popf)


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