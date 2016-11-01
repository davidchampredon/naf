### 
###   ANALYZE SIMULATIONS ALREADY RUN
###

t00 <- as.numeric(Sys.time())

### ==== Load simulation results ====

print('Loading simulation results ...')
load('mc-simul.RData')
t01 <- as.numeric(Sys.time()) ; msgt <- round((t01-t00)/60,2)
print(paste('... simulation results loaded in',msgt,'minutes.'))

source('analysis_tools.R')  
save.plot.to.file <- TRUE


### ==== Merge all MC iterations ====

n.cpu <- parallel::detectCores() - 1

# Merging populations of all MC iterations
# can be very slow, and they very similar
# across MC, so just display the first ('select.mc=1')

if(simul.prm[['baseline_only']]) res.select <- res.list.0
if(!simul.prm[['baseline_only']]) res.select <- res.list

pop   <- merge.pop.mc(res.select,
					  n.cpu = n.cpu, 
					  doparallel = TRUE, 
					  select.mc = 1)

# Merging time series if faster and more informative,
# so it is done across all MC iterations:
ts    <- merge.ts.mc(res.list.0, n.cpu = n.cpu)
tsc   <- merge.ts.mc(res.list.0, n.cpu = n.cpu, is.contact = TRUE)
tssp  <- merge.ts.mc(res.list.0, n.cpu = n.cpu, is.sp = TRUE)

if(exists('res.list')){
	ts.intrv    <- merge.ts.mc(res.list, n.cpu = n.cpu)
	tssp.intrv  <- merge.ts.mc(res.list, n.cpu = n.cpu, is.sp = TRUE)
}

### ==== Plots ====

# World (all social places):
print(' -> Ploting world ...')
if (save.plot.to.file) pdf('plot_world.pdf', width = 10, height = 8)
try( plot.world(res.list.0),  silent = T)
if (save.plot.to.file) dev.off()

# Population:
print(' -> Ploting population ...')
if (save.plot.to.file) pdf('plot_pop.pdf', width = 30, height = 18)
try( plot.population(pop, split.mc=T),  silent = T)
try( plot.n.contacts(tsc),  silent = T)
try( plot.age.contact.matrix.avg(res.list),  silent = T)
try( plot.sp.sz.distrib.new(pop,world.prm) , silent = T)
try( plot.share.same.hh(pop), silent = T)
if (save.plot.to.file) dev.off()

# Time series:
print(' -> Ploting time series ...')
if (save.plot.to.file) pdf('plot_ts.pdf', width = 25,height = 15)
try( plot.epi.timeseries(ts),  silent = T)
try( plot.ts.sp(tssp),  silent = T)
if(exists('res.list')){
	try( plot.epi.timeseries(ts.intrv),  silent = T)
	try( plot.ts.sp(tssp.intrv),  silent = T)
}
if (save.plot.to.file) dev.off()

t1 <- as.numeric(Sys.time())
print(paste('Analysis completed in',round((t1-t00)/60,1),'min'))

