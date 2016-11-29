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
dir.results <- '../results/'

# If outputs for all social places are produced separately:
detailed.analysis <- FALSE

### ==== Merge all MC iterations ====

n.cpu <- parallel::detectCores() - 1

if(simul.prm[['baseline_only']]) res.select <- res.list.0
if(!simul.prm[['baseline_only']]) res.select <- res.list


calc.fizzles <- TRUE

if(calc.fizzles){
	pop.all.mc   <- merge.pop.mc(res.select,
								 n.cpu = n.cpu, 
								 doparallel = TRUE)
}

# Merging populations of all MC iterations
# can be very slow, and they very similar
# across MC, so just display the first non-fizzled('select.mc=idx.mc.no.fizz[1]')	

idx.mc         <- unique(pop.all.mc$mc)
fizz.mc        <- identify.fizzle(pop.all.mc)
fizz.mc.idx    <- which(fizz.mc==TRUE)
idx.mc.no.fizz <- idx.mc[-fizz.mc.idx]

pop <- subset(pop.all.mc, mc==idx.mc.no.fizz[1])

# Merging time series if faster and more informative,
# so it is done across all MC iterations:
res.no.fizz <- list()
for(i in seq_along(idx.mc.no.fizz)) res.no.fizz[[i]] <- res.list.0[[i]]

ts    <- merge.ts.mc(res.no.fizz, n.cpu = n.cpu)
tsc   <- merge.ts.mc(res.no.fizz, n.cpu = n.cpu, is.contact = TRUE)

if(detailed.analysis) tssp  <- merge.ts.mc(res.no.fizz, n.cpu = n.cpu, is.sp = TRUE)

if(exists('res.list')){
	ts.intrv    <- merge.ts.mc(res.list, n.cpu = n.cpu)
	if(detailed.analysis) tssp.intrv  <- merge.ts.mc(res.list, n.cpu = n.cpu, is.sp = TRUE)
}

### ==== Plots ====

# filenames
fname.world <- paste0(dir.results,'plot_world.pdf')
fname.pop   <- paste0(dir.results,'plot_pop.pdf')
fname.ts    <- paste0(dir.results,'plot_ts.pdf')

# World (all social places):
print(' -> Ploting world ...')
if (save.plot.to.file) pdf(fname.world, width = 10, height = 8)
try( plot.world(res.list.0),  silent = T)
if (save.plot.to.file) dev.off()


# Population:
print(' -> Ploting population ...')
if (save.plot.to.file) pdf(fname.pop, width = 30, height = 18)

try( plot.population(pop, split.mc=F),  silent = T)
try( plot.n.contacts(tsc),  silent = T)
try( plot.age.contact.matrix.avg(res.no.fizz),  silent = T)
try( plot.sp.sz.distrib.new(pop,world.prm) , silent = T)
try( plot.share.same.hh(pop), silent = T)
if(calc.fizzles) try( plot.prop.fizzles(pop.all.mc), silent = T)
	
if (save.plot.to.file) dev.off()

# Time series:
print(' -> Ploting time series ...')
if (save.plot.to.file) pdf(fname.ts, width = 25,height = 15)
try( plot.epi.timeseries(ts),  silent = T)
try( plot.ts.sp(tssp),  silent = T)
if(exists('res.list')){
	try( plot.epi.timeseries(ts.intrv),  silent = T)
	try( plot.ts.sp(tssp.intrv),  silent = T)
}
if (save.plot.to.file) dev.off()

t1 <- as.numeric(Sys.time())
print(paste('Analysis completed in',round((t1-t00)/60,1),'min'))

