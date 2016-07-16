##################################################################
######
######    MINIMAL TEST FOR 'stiagent' LIBRARY
######
######
##################################################################


library(naf,lib.loc = "./lib")

t0 <- Sys.time()

prm <- list()

prm[['a']] <- 123

res <- naf_test(prm)

plot(res[[1]],ylim=c(0,100), typ='o')

t1 <- Sys.time()
message(paste("time elapsed:",round(t1-t0,1),"sec"))

