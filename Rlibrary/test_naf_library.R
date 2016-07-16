##################################################################
######
######    MINIMAL TEST FOR 'stiagent' LIBRARY
######
######
##################################################################

library(naf,lib.loc = "./lib")

t0 <- Sys.time()



prm <- list()
prm[['rnd_seed']] <- 123

set.seed(prm[['rnd_seed']])
 
res <- naf_test(prm)

plot(res[[1]],ylim=c(0,30), typ='o')

 t1 <- Sys.time()
 message(paste("time elapsed:",round(t1-t0,1),"sec"))

