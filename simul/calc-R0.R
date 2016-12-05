### 
###   CALCULATE R0 OF SIMULATIONS ALREADY RUN
###

t00 <- as.numeric(Sys.time())

### ==== Load simulation results ====

print('Loading simulation results ...')
load('mc-simul.RData')
t01 <- as.numeric(Sys.time()) ; msgt <- round((t01-t00)/60,2)
print(paste('... simulation results loaded in',msgt,'minutes.'))

source('analysis_tools.R')  

### ==== Merge all MC iterations ====

n.cpu <- parallel::detectCores() - 1
if(simul.prm[['baseline_only']])  res.select <- res.list.0
if(!simul.prm[['baseline_only']]) res.select <- res.list

pop.all.mc   <- merge.pop.mc(res.select,
							 n.cpu = n.cpu, 
							 doparallel = TRUE)

R0 <- calc.R0(pop.all.mc)


