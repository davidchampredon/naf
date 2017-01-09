### 
###   COMPARE SIMULATIONS (BASELINE vs INTERVENTION) ALREADY RUN
###

library(plyr)
library(ggplot2)
theme_set(theme_bw())

compare.simul.scen <- function(scen.id, 
							   dir.save.rdata,
							   dir.results,
							   do.secondary = FALSE) {
	ct0 <- as.numeric(Sys.time())
	library(tidyr)
	library(dplyr)
	
	### ==== Load simulation results ====
	
	print(paste('Loading simulation results for scenario',scen.id,'...'))
	load(paste0(dir.save.rdata,'mc-simul-',scen.id,'.RData'))
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
		if(mx > nrow(z)) {
			tmx <- z0$time 
			tmn <- z$time
		}
		if(mx <= nrow(z)) {
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
					 tot.death = sum(1-is_alive),
					 tot.treat = sum(is_treated)))
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
					tot.treat.baseline = a0$tot.treat,
					# Raw diffences with baseline:
					d.inf   = a$tot.inf   - a0$tot.inf,
					d.sympt = a$tot.sympt - a0$tot.sympt,
					d.hosp  = a$tot.hosp  - a0$tot.hosp,
					d.death = a$tot.death - a0$tot.death,
					d.treat = a$tot.treat - a0$tot.treat)
	
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
	
	# Filter out fizzles, if any:
	threshold.fizzle <- 0.01
	idx.not.fizzle <- which(d$tot.inf.baseline > d$popsize * threshold.fizzle)
	if(nrow(d) > length(idx.not.fizzle))	d <- d[idx.not.fizzle,]
	
	# Save main results:
	result.scen <- d
	save(list = 'result.scen',
		 file = paste0(dir.save.rdata,'result-scen-',scen.id,'.RData'), 
		 compress = FALSE)
	
	
	# ==== Secondary results ====
	if(do.secondary){
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
			 file = paste0(dir.save.rdata,'result-secondary-',scen.id,'.RData'),
			 compress = FALSE)
	}
	# ----- Plots time series -----
	
	pdf(paste0(dir.results,'plot-compare-',scen.id,'.pdf'),
		width = 15, height = 10)
	
	plot.epi.timeseries.comp(u)
	plot.ts.comp.all(dfall)
	
	dev.off()
	t2 <- as.numeric(Sys.time())
	print(paste('Full comparison completed in',
				round((t2-ct0)/60,2),'minutes for scenario',scen.id,'.'))
}

merge.result.scen <- function(scen.id, dir.save.rdata, do.secondary=FALSE) 
{
	# Main results:
	
	bfirst <- TRUE
	for(i in seq_along(scen.id)){
		load(paste0(dir.save.rdata,'result-scen-',scen.id[i],'.RData'))
		if(bfirst)  result.scen.all <- result.scen
		if(!bfirst) result.scen.all <- rbind(result.scen.all, result.scen)
		bfirst <- FALSE
	}
	save(list = 'result.scen.all', 
		 file = paste0(dir.save.rdata,'result-scen-all.RData'),
		 compress = FALSE)
	
	# Secondary results:
	
	secondary.all <- NULL
	if(do.secondary){
		bfirst <- TRUE
		for(i in seq_along(scen.id)){
			load(paste0(dir.save.rdata,'result-secondary-',scen.id[i],'.RData'))
			if(bfirst)  secondary.all <- result.secondary
			if(!bfirst) secondary.all <- rbind(secondary.all, result.secondary)
			bfirst <- FALSE
		}
		save(list = 'secondary.all', 
			 file = paste0(dir.save.rdata,'secondary-all.RData'),
			 compress = FALSE)
	}
	return(list(main      = result.scen.all,
				secondary = secondary.all) )
}

plot.multi.scen.res <- function(result.scen.all, 
								dir,
								file.scen.prm.list){
	idx <- which( substr(names(result.scen.all),1,6)=='rel.d.'  )
	x <- gather(result.scen.all,'type','rel.diff',idx)
	
	pdf(paste0(dir,'plot-compare-all.pdf'), width=20, height=15)
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

plot.rate.reduc <- function(result.scen.all, 
							dir, 
							file.scen.prm.list){
	
	# Retrieve scenarios definitions:
	spl <- read.csv(file.scen.prm.list) # DEBUG::  file.scen.prm.list <- 'scenario-prm-list.csv'
	names(result.scen.all)[names(result.scen.all)=='scenario'] <- 'scenario_id' 
	
	# Merge results and scenarios definition:
	x <- join(result.scen.all, spl, by='scenario_id')
	
	# Use only what's necessary for plot:
	output.stats <- function(x, resp.var) {
		idx <- which(is.nan(x[[resp.var]]) | is.infinite(x[[resp.var]]))
		if(length(idx)>0) x <- x[-idx, ]
		
		z <- ddply(x,
				   c('interv_cvg_rate', 'interv_target', 'contact_rate_mean',
				     'interv_cvg_max_prop','interv_start'),
				   function(df){c(
				   	mn  = mean(-df[[resp.var]]),
				   	md  = quantile(-df[[resp.var]], probs = 0.50,names = F, na.rm = T),
				   	qlo = quantile(-df[[resp.var]], probs = 0.05,names = F, na.rm = T),
				   	qhi = quantile(-df[[resp.var]], probs = 0.95,names = F, na.rm = T),
				   	n   = length(df[[resp.var]]))
				   }
		)
		return(z)
	}
	
	z.inf   <- output.stats(x, resp.var = 'rel.d.inf')
	z.sympt <- output.stats(x, resp.var = 'rel.d.sympt')
	z.treat <- output.stats(x, resp.var = 'rel.d.treat')
	z.hosp  <- output.stats(x, resp.var = 'rel.d.hosp')
	z.death <- output.stats(x, resp.var = 'rel.d.death')
	
	plot.reduc.curve <- function(z, title='')
	{
		z$key <- paste0(z$interv_target,
						' ; CR=',z$contact_rate_mean,
						' ; CvgMx=',z$interv_cvg_max_prop,
						' ; Start= ',z$interv_start)
		
		g <- ggplot(z,aes(x=interv_cvg_rate, y=mn, colour = factor(interv_start)), alpha=0.6) 
		g <- g + geom_point(size=1) + geom_line(size=0.5) + geom_linerange(aes(ymin=qlo, ymax=qhi), size=5, alpha=0.2)
		g <- g + geom_point(data=z, aes(x=interv_cvg_rate, y=md, colour = factor(interv_start)), size=1, shape=3)
		g <- g + facet_wrap(~key)
		g <- g + scale_color_brewer(palette = 'Dark2') + ylab('Mean relative reduction') 
		g <- g + ggtitle(title)
		plot(g)
	}
	
	plot.start.compare <- function(z, title='') {
		z$key <- paste0(z$interv_target,
						' ; CR=',z$contact_rate_mean,
						' ; CvgMx=',z$interv_cvg_max_prop)
		
		
		g <- ggplot(z) + geom_line(aes(x=interv_cvg_rate, y=mn, colour=factor(interv_start)), size=3, alpha=0.7)
		g <- g + geom_point(aes(x=interv_cvg_rate, y=mn, colour=factor(interv_start)), size=4, alpha=0.7)
		g <- g + facet_wrap(~key)
		g <- g + scale_color_brewer(palette = 'BrBG')
		g <- g + ggtitle(title) + ylab('Mean relative reduction')
		plot(g)
	}
	
	plot.mc.cr <- function(x){
		fiz <- ddply(x,c('scenario_id'),summarize, n=length(mc), cr = mean(contact_rate_mean))
		g <- ggplot(fiz,aes(x=factor(cr), y=n)) + geom_boxplot()
		g <- g + ggtitle('Number of MC iterations') + xlab('Mean contact rate')
		plot(g)
	}
	
	
	
	
	pdf(paste0(dir,'plot-rate-reduc.pdf'), width=18, height=17)
	plot.mc.cr(x)
	
	zlist <- list(z.inf,z.sympt,z.hosp,z.death)
	titlelist <- list('All Infections','Symptomatic Infections','Hospitalized','Deaths')
	
	for(i in seq_along(zlist)) plot.reduc.curve(zlist[[i]], title = titlelist[[i]])
	for(i in seq_along(zlist)) plot.start.compare(zlist[[i]], title = titlelist[[i]])
	
	dev.off()
}





