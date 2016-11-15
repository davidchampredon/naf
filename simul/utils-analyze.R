

analyze.simul.scen <- function(scen.id,	
							   dir.save.rdata,
							   dir.results,
							   save.plot.to.file = TRUE,
							   detailed.analysis = FALSE) {
	
	t00 <- as.numeric(Sys.time())
	
	# Load simulation results
	print(paste('Loading simulation results for scenario',scen.id,'...'))
	load(paste0(dir.save.rdata,'mc-simul-',scen.id, '.RData'))
	t01 <- as.numeric(Sys.time()) ; msgt <- round((t01-t00)/60,2)
	print(paste('... simulation results for scenario',scen.id,'loaded in',msgt,'minutes.'))
	
	source('analysis_tools.R')  
	
	### ==== Merge all MC iterations ====
	
	n.cpu <- parallel::detectCores() - 1
	
	# Merging populations of all MC iterations
	# can be very slow, and they very similar
	# across MC, so just display the first ('select.mc=1')
	pop   <- merge.pop.mc(res.list.0,
						  n.cpu = n.cpu, 
						  doparallel = TRUE, 
						  select.mc = 1)
	
	# Merging time series if faster and more informative,
	# so it is done across all MC iterations:
	ts    <- merge.ts.mc(res.list.0, n.cpu = n.cpu)
	tsc   <- merge.ts.mc(res.list.0, n.cpu = n.cpu, is.contact = TRUE)
	if(detailed.analysis) tssp  <- merge.ts.mc(res.list.0, n.cpu = n.cpu, is.sp = TRUE)
	
	if(exists('res.list')){
		ts.intrv    <- merge.ts.mc(res.list, n.cpu = n.cpu)
		if(detailed.analysis) tssp.intrv  <- merge.ts.mc(res.list, n.cpu = n.cpu, is.sp = TRUE)
	}
	
	### ==== Plots ====
	
	# World (all social places):
	print(' -> Ploting world ...')
	if (save.plot.to.file) pdf(paste0(dir.results,'plot_world_',scen.id,'.pdf'), 
							   width = 10, height = 8)
	try( plot.world(res.list.0),  silent = T)
	if (save.plot.to.file) dev.off()
	
	# Population:
	print(' -> Ploting population ...')
	if (save.plot.to.file) pdf(paste0(dir.results,'plot_pop_',scen.id,'.pdf'), 
							   width = 30, height = 18)
	try( plot.population(pop, split.mc=T),  silent = T)
	try( plot.n.contacts(tsc),  silent = T)
	try( plot.age.contact.matrix.avg(res.list),  silent = T)
	try( plot.sp.sz.distrib.new(pop,world.prm) , silent = T)
	try( plot.share.same.hh(pop), silent = T)
	if (save.plot.to.file) dev.off()
	
	# Time series:
	print(' -> Ploting time series ...')
	if (save.plot.to.file) pdf(paste0(dir.results,'plot_ts_',scen.id,'.pdf'), 
							   width = 25,height = 15)
	try( plot.epi.timeseries(ts),  silent = T)
	if(detailed.analysis) try( plot.ts.sp(tssp),  silent = T)
	if(exists('res.list')){
		try( plot.epi.timeseries(ts.intrv),  silent = T)
		if(detailed.analysis) try( plot.ts.sp(tssp.intrv),  silent = T)
	}
	if (save.plot.to.file) dev.off()
	
	t1 <- as.numeric(Sys.time())
	msg.a <- paste('Analysis completed in',round((t1-t00)/60,1),'min')
	print(msg.a)
	message(msg.a)
}