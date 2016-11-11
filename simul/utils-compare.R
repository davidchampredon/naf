### 
###   COMPARE SIMULATIONS (BASELINE vs INTERVENTION) ALREADY RUN
###

library(plyr)
library(ggplot2)


compare.simul.scen <- function(scen.id) {
	
	ct0 <- as.numeric(Sys.time())
	library(tidyr)
	library(dplyr)
	
	### ==== Load simulation results ====
	
	print(paste('Loading simulation results for scenario',scen.id,'...'))
	load(paste0('mc-simul-',scen.id,'.RData'))
	t1 <- as.numeric(Sys.time())
	print(paste('... simulation results for scenario',scen.id,'loaded in',round((t1-ct0)/60,1),'minutes.'))
	
	n.mc  <- length(res.list)
	print(paste('Number of MC iterations:',n.mc))
	
	source('analysis_tools.R')  
	save.plot.to.file <- TRUE
	max.cpu <- 2
	
	
	# Merge all MC iterations
	
	n.cpu <- min(parallel::detectCores(), max.cpu)
	ts    <- merge.ts.mc(res.list,   n.cpu = n.cpu)
	ts0   <- merge.ts.mc(res.list.0, n.cpu = n.cpu)
	
	ts0$scen <- 'baseline'
	ts$scen  <- 'interv'
	u        <- rbind.data.frame(ts0,ts)
	
	tmp     <- list()
	
	# Loop to calculate differences
	# b/w baseline and intervention:
	for(i in 1:n.mc){
		print(paste('Calculating scenario differences',i,'/',n.mc))
		z0 <- subset(ts0,mc==i)
		z  <- subset(ts, mc==i)
		
		pop.size <- z$nS[1]+z$nE[1]+z$nIa[1]+z$nIs[1]
		
		mx <- max(nrow(z),nrow(z0))
		mn <- min(nrow(z),nrow(z0))
		if(mx>nrow(z)) {
			tmx <- z0$time 
			tmn <- z$time
		}
		else {
			tmx <- z$time
			tmn <- z0$time
		}
		
		tmp[[i]] <- data.frame(time = tmx, 
							   mc   = factor(i),
							   dinc = diff.flex(z$incidence, z0$incidence)/pop.size ,
							   dIa  = diff.flex(z$nIa, z0$nIa)/pop.size,
							   dIs  = diff.flex(z$nIs, z0$nIs)/pop.size
		)
	}
	dfall <- do.call('rbind.data.frame',tmp)
	
	# ----- Comparison tables -----
	
	# ==== Main results ====
	pop0 <- merge.pop.mc(res.list = res.list.0, 
						 n.cpu = 2,
						 doparallel = TRUE,
						 select.mc = 1:n.mc)
	pop <- merge.pop.mc(res.list = res.list, 
						n.cpu = 2,
						doparallel = TRUE,
						select.mc = 1:n.mc)
	
	main.results <- function(df){
		return(ddply(df,c('mc'),summarise, 
					 tot.pop   = length(unique(id_indiv)),
					 tot.inf   = sum(is_recovered),
					 tot.sympt = sum(was_symptomatic),
					 tot.hosp  = sum(was_hosp),
					 tot.death = sum(1-is_alive)))
	}
	a0 <- main.results(pop0)
	a0$intervention <- 'baseline'
	a <- main.results(pop)
	a$intervention <- 'interv' #paste('interv',scen.id,sep='-')
	
	# Data frame of raw differences:
	d <- data.frame(scenario = scen.id,
					mc = a$mc,
					popsize   = a0$tot.pop,
					tot.inf.baseline   = a0$tot.inf,
					tot.sympt.baseline = a0$tot.sympt,
					tot.hosp.baseline  = a0$tot.hosp,
					tot.death.baseline = a0$tot.death,
					# Raw diffences with baseline:
					d.inf   = a$tot.inf   - a0$tot.inf,
					d.sympt = a$tot.sympt - a0$tot.sympt,
					d.hosp  = a$tot.hosp  - a0$tot.hosp,
					d.death = a$tot.death - a0$tot.death)
	
	# Add _relative_ differences calculations
	n <- ncol(d)
	idx.diff     <- which( substr(names(d),1,2)=='d.'  )
	idx.tot.base <- which( substr(names(d),1,4)=='tot.' )
	stopifnot(length(idx.diff)==length(idx.tot.base))
	
	for(j in 1:length(idx.diff)) {
		v <- d[,idx.diff[j] ] / d[,idx.tot.base[j]]
		d <- cbind(d, v)
		nam <- paste0('rel.',names(d)[ idx.diff[j] ])
		names(d)[ncol(d)] <- nam
	}
	# Save main results:
	result.scen <- d
	save(list = 'result.scen',
		 file = paste0('result-scen-',scen.id,'.RData'), 
		 compress = FALSE)
	
	
	# ==== Secondary results ====
	
	secondary.results <- function(df){
		select(df,age,
			   immunity_hum, immunity_cell, 
			   frailty,
			   doi, n_secondary_cases)
	}
	b0 <- secondary.results(pop0)
	b0$intervention <- 'baseline'
	b0$scenario <- scen.id
	b <- secondary.results(pop)
	b$intervention <- 'interv' #paste('interv',scen.id,sep='-')
	b$scenario <- scen.id
	
	# Save secondary results
	result.secondary <- rbind(b0,b)
	save(list = 'result.secondary',
		 file = paste0('result-secondary-',scen.id,'.RData'),
		 compress = FALSE)
	
	
	# ----- Plots time series -----
	
	pdf(paste0('plot-compare-',scen.id,'.pdf'),
		width = 15, height = 10)
	
	plot.epi.timeseries.comp(u)
	plot.ts.comp.all(dfall)
	
	dev.off()
	t2 <- as.numeric(Sys.time())
	print(paste('Full comparison completed in',
				round((t2-ct0)/60,2),'minutes for scenario',scen.id,'.'))
}

merge.result.scen <- function(scen.id) {
	
	# Main results:
	for(i in seq_along(scen.id)){
		load(paste0('result-scen-',scen.id[i],'.RData'))
		if(i==scen.id[1]) result.scen.all <- result.scen
		else result.scen.all <- rbind(result.scen.all, result.scen)
	}
	save(list = 'result.scen.all', 
		 file = 'result-scen-all.RData',
		 compress = FALSE)
	
	# Secondary results:
	for(i in seq_along(scen.id)){
		load(paste0('result-secondary-',scen.id[i],'.RData'))
		if(i==scen.id[1]) secondary.all <- result.secondary
		else secondary.all <- rbind(secondary.all, result.secondary)
	}
	save(list = 'secondary.all', 
		 file = 'secondary-all.RData',
		 compress = FALSE)
	
	return(list(main = result.scen.all,
				secondary = secondary.all) )
}

plot.multi.scen.res <- function(result.scen.all, 
								dir,
								file.scen.prm.list){
	idx <- which( substr(names(result.scen.all),1,6)=='rel.d.'  )
	x <- gather(result.scen.all,'type','rel.diff',idx)
	
	pdf(paste0(dir,'plot-compare-all.pdf'), width=12, height=8)
	g <- ggplot(x) + geom_boxplot(aes(x=factor(type),y=rel.diff, 
									  fill=factor(scenario)))
	g <- g + geom_hline(yintercept=0, colour='black',linetype=2)
	plot(g)
	
	g1 <- g + facet_wrap(~scenario)
	plot(g1)
	
	g2 <- ggplot(x) + geom_boxplot(aes(x=factor(scenario),y=rel.diff, 
									   fill=factor(scenario)))+ facet_wrap(~type, scales = 'free_y')
	g2 <- g2 + geom_hline(yintercept=0, colour='black',linetype=2)
	plot(g2)
	
	# - - - - Change this plot - - - - 
	# Retrieve scenarios definitions:
	spl <- read.csv(file.scen.prm.list) # DEBUG::  file.scen.prm.list <- 'scenario-prm-list.csv'
	names(result.scen.all)[names(result.scen.all)=='scenario'] <- 'scenario_id' 
	# Merge results and scenarios definition:
	x <- join(result.scen.all, spl, by='scenario_id')
	# Plots
	g3 <- ggplot(x)
	g3 <- g3 + geom_point(aes(x = interv_cvg_rate,
							  y = -rel.d.sympt,
							  colour = factor(interv_start)), alpha=0.6)
	g3 <- g3 + facet_wrap(~interv_target+contact_rate_mean+interv_cvg_max_prop)
	plot(g3)
	
	dev.off()
}

plot.secondary.res <- function(df, dir){
	# df <- res.all[['secondary']]
	# df$scenario <- df$scen
	age.bucket <- 2
	df$age.group <- round(df$age/age.bucket)*age.bucket
	
	x <- ddply(df, c('scenario','age.group','intervention'),
			   summarize, 
			   immh = mean(immunity_hum),
			   immc = mean(immunity_cell),
			   frailty = mean(frailty))
	
	g.age <- ggplot(df) + geom_density(aes(x=age.group,
										   colour=factor(scenario),
										   linetype = factor(intervention)))
	
	g.immh <- ggplot(x) + geom_line(aes(x=age.group,y=immh,
								   colour=factor(scenario),
								   linetype = factor(intervention)))
	
	g.immc <- ggplot(x) + geom_line(aes(x=age.group,y=immc,
									colour=factor(scenario),
									linetype = factor(intervention)))
	
	g.frail <- ggplot(x) + geom_line(aes(x=age.group, y=frailty,
										 colour=factor(scenario),
										 linetype = factor(intervention)))
	
	
	df.transm <- subset(df,n_secondary_cases>0)
	a <- ddply(df.transm,c('scenario','intervention'),summarise,
			   R.m   = mean(n_secondary_cases),
			   doi.m = mean(doi))
	a$scenario <- as.factor(a$scenario)
	a$intervention <- as.factor(a$intervention)
	
	g.R <- ggplot(a) + geom_point(aes(y = R.m,
									  x = scenario,
									  shape = intervention))
	g.R <- g.R + ggtitle('Mean number secondary cases')
	
	g.doi <- ggplot(a) + geom_point(aes(y = doi.m,
										x = scenario,
										shape = intervention))
	g.doi <- g.doi + ggtitle('Mean DOI')
	
	pdf(paste0(dir,'plot-secondary.pdf'), width=15, height = 10)
	grid.arrange(g.immh,g.immc, g.frail)
	grid.arrange(g.R, g.doi)
	grid.arrange(g.age)
	dev.off()
}

plot2 <- function(result.scen.all, dir, file.scen.prm.list){
	
	# Retrieve scenarios definitions:
	spl <- read.csv(file.scen.prm.list) # DEBUG::  file.scen.prm.list <- 'scenario-prm-list.csv'
	names(result.scen.all)[names(result.scen.all)=='scenario'] <- 'scenario_id' 
	
	# Merge results and scenarios definition:
	x <- join(result.scen.all, spl, by='scenario_id')
	
	# Plots
	g <- ggplot(x)
	
	g <- g + geom_point(aes(x = interv_cvg_rate,
								   y = -rel.d.sympt,
								   colour = factor(interv_start)), alpha=0.6)
	
	# g <- g + geom_boxplot(aes(x=factor(interv_cvg_rate), 
	# 						  y = -rel.d.sympt,
	# 						  colour = factor(interv_start)))
	
	g <- g + facet_wrap(~interv_target+contact_rate_mean+interv_cvg_max_prop)
	plot(g)
}