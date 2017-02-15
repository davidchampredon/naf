### 
###   COMPARE SIMULATIONS (BASELINE vs INTERVENTION) ALREADY RUN
###

library(plyr)
library(ggplot2) ; theme_set(theme_bw())
library(tidyr)
library(dplyr)

#' Define age groups
ageGroup <- function(x) {
    age.1 <- 5
    age.2 <- 18
    age.3 <- 65
    x$ageGroup <- NA
    x$ageGroup[ x$age < age.1] <- paste0('0_',age.1)
    x$ageGroup[ age.1 <= x$age & x$age < age.2] <- paste(age.1,age.2, sep='_')
    x$ageGroup[ age.2 <= x$age & x$age < age.3] <- paste(age.2,age.3, sep='_')
    x$ageGroup[ age.3 <= x$age ] <- paste0(age.3,'_over')
    return(x)
}

#' Summarize main outcomes
main.results <- function(df, intervention.type, do.ageGroup = FALSE){
    
    sel <- 'mc'
    if(do.ageGroup) sel <- c('mc','ageGroup')
    
    res <- ddply(df,sel,summarise, 
                 tot.pop   = length(unique(id_indiv)),
                 tot.inf   = sum(is_recovered),
                 tot.sympt = sum(was_symptomatic),
                 tot.hosp  = sum(was_hosp),
                 tot.death = sum(1-is_alive),
                 tot.treat = sum(is_treated))
    res$intervention <- intervention.type
    return(res)
}

#' Calculate raw difference of main outcomes with baseline
diff_df <- function(a0, a, scen.id, do.ageGroup = FALSE) {
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
    if(do.ageGroup) d$ageGroup <- a$ageGroup
    return(d)
}

#' Calculate relative difference of main outcomes with baseline
reldiff_df <- function(d) {
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
    return(d)
}


#' Filter out fizzles, if any:
filter.fizzle <- function(d) {
    threshold.fizzle <- 0.01
    idx.not.fizzle <- which(d$tot.inf.baseline > d$popsize * threshold.fizzle)
    if(nrow(d) > length(idx.not.fizzle))
        d <- d[idx.not.fizzle,]
    return(d)
}

compare.simul.scen <- function(scen.id, 
							   dir.save.rdata,
							   dir.results,
							   do.secondary = FALSE) {
	ct0 <- as.numeric(Sys.time())

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
	
	pop0 <- merge.pop.mc(res.list = res.list.0, 
						 n.cpu = 2,
						 doparallel = FALSE,
						 select.mc = 1:n.mc)
	pop <- merge.pop.mc(res.list = res.list, 
						n.cpu = 2,
						doparallel = FALSE,
						select.mc = 1:n.mc)

	# Assign age groups:
	pop0 <- ageGroup(pop0)
	pop  <- ageGroup(pop)
	
	# Summarize the main outcomes:
	a0 <- main.results(df= pop0,'baseline')
	a  <- main.results(df= pop,'interv')
	a0.ageGroup <- main.results(df= pop0, 'baseline', do.ageGroup = TRUE)
	a.ageGroup  <- main.results(df= pop, 'interv', do.ageGroup = TRUE)

	# Raw difference with baseline:
	d          <- diff_df(a0, a, scen.id, do.ageGroup = FALSE)
	d.ageGroup <- diff_df(a0.ageGroup, a.ageGroup, scen.id, do.ageGroup = TRUE)
	
	# Add _relative_ differences calculations
	d <- reldiff_df(d)
	d.ageGroup <- reldiff_df(d.ageGroup)
	
	# Filter fizzles:
	d <- filter.fizzle(d)
	d.ageGroup <- filter.fizzle(d.ageGroup)
	
	# Save main results:
	result.scen <- d
	result.scen.ageGroup <- d.ageGroup
	
	save(list = c('result.scen', 'result.scen.ageGroup'),
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
		if(bfirst)  {
		    result.scen.all <- result.scen
		    result.scen.all.ageGroup <- result.scen.ageGroup
		}
		if(!bfirst) {
		    result.scen.all <- rbind(result.scen.all, result.scen)
		    result.scen.all.ageGroup <- rbind(result.scen.all.ageGroup, result.scen.ageGroup)
		}
		bfirst <- FALSE
	}
	save(list = c('result.scen.all', 'result.scen.all.ageGroup'), 
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

explicit_variable <- function(var, df){
    df[[var]] <- paste(var,'=',df[[var]])
    return(df)
}

explicit_variable_2 <- function(var, df, newvar){
    df[[newvar]] <- paste(newvar,'=',df[[var]])
    return(df)
}



explicit_all <- function(z){
    
    var.vec <- c('interv_target', 'contact_rate_mean','interv_cvg_max_prop',
                 'interv_efficacy','imm_hum_baseline')
    for(i in seq_along(var.vec)){
        z <- explicit_variable(var = var.vec[i], z)
    }
    return(z)
}

explicit_all_2 <- function(z){
    
    idx <- which(names(z) %in% c('mn', 'md', 'qlo', 'qhi', 'n'))
    var.vec <- c(1:length(names(z)))[-idx]
    
    for(i in seq_along(var.vec)){
        z <- explicit_variable(var = names(z)[var.vec[i]], z)
    }
    return(z)
}


# Use only what's necessary for plot:
output.stats <- function(x, resp.var, do.ageGroup = FALSE) {
    idx <- which(is.nan(x[[resp.var]]) | is.infinite(x[[resp.var]]))
    if(length(idx)>0) x <- x[-idx, ]
    
    sel <- c('interv_cvg_rate', 'interv_target', 'contact_rate_mean',
             'interv_cvg_max_prop','interv_start', 'interv_efficacy',
             'imm_hum_baseline')
    
    if(do.ageGroup) sel <- c(sel,'ageGroup')
    
    z <- ddply(x, sel,
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

plot.reduc.curve <- function(z, title='', do.ageGroup = FALSE)
{
    g <- ggplot(z,aes(x=interv_cvg_rate, y=mn, colour = factor(interv_start)), alpha=0.3) 
    g <- g + geom_point(size=1) + 
        geom_line(size=0.5) + 
        geom_linerange(aes(ymin=qlo, ymax=qhi), size=5, alpha=0.2)
    
    if(!do.ageGroup) g <- g + facet_grid(contact_rate_mean ~ interv_target+interv_efficacy+interv_cvg_max_prop+imm_hum_baseline)
    if(do.ageGroup)  g <- g + facet_grid(contact_rate_mean + ageGroup ~ interv_target+interv_efficacy+interv_cvg_max_prop+imm_hum_baseline)
    g <- g + ggtitle(title) + 
        scale_color_brewer(palette = 'Dark2') + 
        ylab('Mean relative reduction') +
        guides(colour=guide_legend(title="Vaccination lag"))
    plot(g)
}

plot.start.compare <- function(z, title='', do.ageGroup = FALSE) {
    g <- ggplot(z) + geom_line(aes(x=interv_cvg_rate, y=mn, colour=factor(interv_start)), 
                               size=3, alpha=0.7)
    g <- g + geom_point(aes(x=interv_cvg_rate, y=mn, colour=factor(interv_start)), 
                        size=4, alpha=0.7)
    
    if(!do.ageGroup) g <- g + facet_grid(interv_target~contact_rate_mean+interv_efficacy+imm_hum_baseline)
    if(do.ageGroup)  g <- g + facet_grid(interv_target + ageGroup ~ contact_rate_mean+interv_efficacy+imm_hum_baseline)
    
    g <- g + scale_color_brewer(palette = 'BrBG')+ 
        ggtitle(title) + ylab('Mean relative reduction')+ 
        guides(colour=guide_legend(title="Vaccination Lag"))
    plot(g)
}


plot.immhum.compare <- function(z, title='', do.ageGroup = FALSE) {
    g <- ggplot(z) + geom_line(aes(x=interv_cvg_rate, y=mn, 
                                   colour=factor(imm_hum_baseline)), 
                               size=3, alpha=0.7)+ 
        geom_point(aes(x=interv_cvg_rate, y=mn, 
                       colour=factor(imm_hum_baseline)), 
                   size=4, alpha=0.7)
    
    if(!do.ageGroup) g <- g + facet_grid(contact_rate_mean+interv_target ~ interv_efficacy+interv_start)
    if(do.ageGroup)  g <- g + facet_grid(contact_rate_mean+interv_target+ageGroup ~ interv_efficacy+interv_start)
    
    g <- g + #scale_color_brewer(palette = 'BrBG')+ 
        ggtitle(title) + ylab('Mean relative reduction')+ 
        guides(colour=guide_legend(title="Baseline Hum. Imm."))
    plot(g)
}


plot.target.compare <- function(z, title='', do.ageGroup = FALSE) {
    
    g <- ggplot(z) + geom_line(aes(x=interv_cvg_rate, y=mn, 
                                   colour=factor(interv_target) ), 
                               size=3, alpha=0.7)
    g <- g + geom_point(aes(x=interv_cvg_rate, y=mn, 
                            colour=factor(interv_target)), 
                        size=4, alpha=0.7)
    if(!do.ageGroup) g <- g + facet_grid(contact_rate_mean+imm_hum_baseline ~ interv_efficacy+interv_start)
    if(do.ageGroup) g <- g + facet_grid(contact_rate_mean+imm_hum_baseline+ageGroup ~ interv_efficacy+interv_start)
    
    g <- g + guides(colour=guide_legend(title="Vaccination strategy"))
    g <- g + ggtitle(title) + ylab('Mean relative reduction')
    plot(g)
}

plot.efficacy.compare <- function(z, title='', do.ageGroup = FALSE) {
    
    g <- ggplot(z) + 
        geom_line(aes(x=interv_cvg_rate, y=mn, 
                      colour=factor(interv_efficacy) ), 
                  size=3, alpha=0.7)+ 
        geom_point(aes(x=interv_cvg_rate, y=mn, 
                       colour=factor(interv_efficacy)), 
                   size=4, alpha=0.7) + 
        guides(colour=guide_legend(title="Efficacy")) + 
        ggtitle(title) + ylab('Mean relative reduction')
    
    if(!do.ageGroup) g <- g + facet_grid(contact_rate_mean+ imm_hum_baseline ~ interv_target + interv_start)
    if(do.ageGroup)  g <- g + facet_grid(contact_rate_mean+ imm_hum_baseline+ ageGroup ~ interv_target + interv_start)
    plot(g)
    
    effvec <- unique(z$interv_efficacy)
    if(length(effvec)>1){
        a1 <- subset(z, interv_efficacy==effvec[1])
        a2 <- subset(z, interv_efficacy==effvec[2])
        a <- a1
        a$reldiff_eff <- a1$mn/a2$mn -1
        a$diff_eff <- a1$mn-a2$mn
        gd <- ggplot(a) +
            geom_point(aes(x=interv_cvg_rate, y=diff_eff), shape='A',size=4) + geom_line(aes(x=interv_cvg_rate, y=diff_eff),alpha=0.3) + 
            geom_point(aes(x=interv_cvg_rate, y=reldiff_eff),shape='R',size=4) + geom_line(aes(x=interv_cvg_rate, y=reldiff_eff),alpha=0.3) + 
            ggtitle(paste(title,'- Impact of vax efficacy')) + ylab('Relative difference of mean reduction')
        
        if(!do.ageGroup) gd <- gd + facet_grid(contact_rate_mean+ imm_hum_baseline ~ interv_target + interv_start)
        if(do.ageGroup)  gd <- gd + facet_grid(contact_rate_mean+ imm_hum_baseline+ ageGroup ~ interv_target + interv_start)
        plot(gd)
    }
}

plot.cvgmax.compare <- function(z, title='', do.ageGroup = FALSE) {
    z$key <- paste0(z$interv_start,
                    ' ; CR=',z$contact_rate_mean,
                    ' ; CvgMx=',z$interv_cvg_max_prop)
    
    g <- ggplot(z) + geom_line(aes(x=interv_cvg_rate, y=mn, colour=factor(interv_cvg_max_prop)), 
                               size=3, alpha=0.2)
    g <- g + geom_point(aes(x=interv_cvg_rate, y=mn, colour=factor(interv_cvg_max_prop)), 
                        size=4, alpha=0.2)
    if(!do.ageGroup) g <- g + facet_grid(contact_rate_mean~interv_start+interv_target)
    if(do.ageGroup)  g <- g + facet_grid(contact_rate_mean+ageGroup~interv_start+interv_target)
    
    g <- g + guides(colour=guide_legend(title="Max coverage"))
    g <- g + ggtitle(title) + ylab('Mean relative reduction')
    plot(g)
}

plot.mc.cr <- function(x){
    mc.max <- max(unique(x$mc))
    fiz <- ddply(x,c('scenario_id'),summarize, n=length(mc), cr = mean(contact_rate_mean))
    g <- ggplot(fiz,aes(x=factor(cr), y=n)) + geom_boxplot(fill='grey')
    g <- g + ggtitle('Number of MC iterations') + xlab('Mean contact rate') + ylab('Non-fizzled MC')
    g <- g + geom_hline(yintercept=mc.max, linetype=2, colour='tomato', size=1)
    g <- g + coord_cartesian(ylim=c(0,mc.max*1.1))
    plot(g)
}

plot.rate.reduc <- function(df, 
                            dir, 
                            file.scen.prm.list,
                            do.ageGroup){
    
    # Retrieve scenarios definitions:
    spl <- read.csv(file.scen.prm.list) # DEBUG::  file.scen.prm.list <- 'scenario-prm-list.csv'
    names(df)[names(df)=='scenario'] <- 'scenario_id' 
    
    # Merge results and scenarios definition:
    x <- join(df, spl, by='scenario_id')
    
    z.inf   <- output.stats(x, resp.var = 'rel.d.inf', do.ageGroup = do.ageGroup)
    z.sympt <- output.stats(x, resp.var = 'rel.d.sympt', do.ageGroup = do.ageGroup)
    z.treat <- output.stats(x, resp.var = 'rel.d.treat', do.ageGroup = do.ageGroup)
    z.hosp  <- output.stats(x, resp.var = 'rel.d.hosp', do.ageGroup = do.ageGroup)
    z.death <- output.stats(x, resp.var = 'rel.d.death', do.ageGroup = do.ageGroup)
    
    z.inf   <- explicit_all(z.inf)
    z.sympt <- explicit_all(z.sympt)
    z.hosp  <- explicit_all(z.hosp)
    z.death <- explicit_all(z.death)
    
    # Plot and save :
    figname <- paste0('plot-rate-reduc',
                      ifelse(do.ageGroup,'-ageGroup',''),
                      '.pdf')
    pdf(paste0(dir,figname), width=30, height=17)
    plot.mc.cr(x)
    
    zlist <- list(z.inf, 
                  z.sympt, 
                  z.hosp, 
                  z.death)
    titlelist <- list('All Infections',
                      'Symptomatic Infections',
                      'Hospitalized',
                      'Deaths')
    # Flag that tests if there were more than one value
    # for a given parameter. 
    # If two or more values, plot comparison.
    u_immh   <- length(unique(spl$imm_hum_baseline))
    u_cvgmax <- length(unique(spl$interv_cvg_max_prop))
    
    for(i in seq_along(zlist)) plot.reduc.curve(zlist[[i]], title = titlelist[[i]], do.ageGroup = do.ageGroup)
    for(i in seq_along(zlist)) plot.start.compare(zlist[[i]], title = titlelist[[i]], do.ageGroup = do.ageGroup)
    for(i in seq_along(zlist)) plot.target.compare(zlist[[i]], title = titlelist[[i]], do.ageGroup = do.ageGroup)
    for(i in seq_along(zlist)) plot.efficacy.compare(zlist[[i]], title = titlelist[[i]], do.ageGroup = do.ageGroup)
    if(u_immh>1) 
        for(i in seq_along(zlist)) plot.immhum.compare(zlist[[i]], title = titlelist[[i]],do.ageGroup = do.ageGroup)
    if(u_cvgmax>1)
        for(i in seq_along(zlist)) plot.cvgmax.compare(zlist[[i]], title = titlelist[[i]],do.ageGroup = do.ageGroup)
    
    dev.off()
}

format.df.for.figure <- function(z){
    
    z$interv_cvg_rate <- z$interv_cvg_rate * 1e5
    
    z$`Vaccination scenario` <- 'Random'
    z$`Vaccination scenario`[z$interv_target=='young_old'] <- 'Young & Senior'
    
    z$`Vaccination lag` <- round(z$interv_start/7)
    
    z$`Vaccine efficacy` <- z$interv_efficacy
    
    if(length(unique(z$interv_cvg_max_prop))==1) 
        z <- subset(z, select= -interv_cvg_max_prop)
    
    # TO DO: CHANGE THAT, AUTOMATE USING MODEL OUTPUTS!
    z$`Final Cum. Incidence` <- 0.5
    z$`Final Cum. Incidence`[z$contact_rate_mean == max(z$contact_rate_mean)] <- 0.7
    
    return(z)
}

figure.TODELETE <- function(z, title='') {
    
    z <- subset(z, interv_efficacy==0.8)
    #z <- subset(z, imm_hum_baseline==0.1)
    z <- subset(z, `Final Cum. Incidence`==0.5)
    
    z <- explicit_variable(var ='Vaccination scenario', df = z )
    z <- explicit_variable(var ='Final Cum. Incidence', df = z )
    
    g <- ggplot(z) 
    
    g <- g + geom_line(aes(x=interv_cvg_rate, 
                           y=mn, 
                           colour=factor(`Vaccination lag`)), 
                       size=3, alpha=0.6)+ 
        geom_errorbar(aes(x = interv_cvg_rate , 
                            ymin=qlo, 
                            ymax=qhi, 
                            colour=factor(`Vaccination lag`)),
                        alpha = 0.3, size=1,
                        #position='dodge',
                        width = 25) +
        geom_point(aes(x=interv_cvg_rate , 
                       y=mn, 
                       colour=factor(`Vaccination lag`)), 
                   size=4, alpha=1)
    
    g <- g + facet_grid(~`Vaccination scenario`)
    
    g <- g + scale_x_continuous(breaks=unique(z$interv_cvg_rate),)
    
    g <- g + scale_color_brewer(palette = 'BrBG')+ 
        ggtitle(title) + ylab('Mean relative reduction')+ 
        xlab('Vaccine administration rate (per 100,000 per day)') + 
        guides(colour=guide_legend(title="Vaccination Lag"))
    
    plot(g)
}


figure.seyed <- function(zlist) {
    
    df <- do.call('rbind.data.frame', zlist)
    
    df$outcome  <- factor(df$outcome, levels = c("Symptomatic Infections", "Hospitalized", "Deaths"))
    df$ageGroup <- factor(df$ageGroup, levels = c("0_5", "5_18", "18_65","65_over"))
    
    cr.vec <- unique(df$contact_rate_mean)
    ve.vec <- unique(df$`Vaccine efficacy`)
    
    q <- expand.grid(cr.vec, ve.vec)
    
    subplot <- function(cr, ve) {
        df <- subset(df, contact_rate_mean==cr & `Vaccine efficacy`==ve)
        
        g <- ggplot(df, aes(x=interv_cvg_rate, y=mn, color=factor(interv_start))) +
            geom_line(alpha=0.5, size=2) +
            geom_point(alpha=0.5, size=4) +
            geom_point(size=4, shape=1) +
            facet_grid(ageGroup ~  outcome) +
            guides(color=guide_legend(title="Vacc. Lag (days)")) +
            coord_cartesian(ylim=c(0,1)) + 
            ggtitle(paste('CR =',cr,'\n VaxEff =',ve))+
            scale_colour_brewer(palette = 'RdYlGn')+
            xlab("Vaccine administration rate (per 100,000 per day)") + ylab("Mean relative reduction") +
            theme_gray()
        plot(g)
    }
    
    subplot.2 <- function(cr, ve) {
        df <- subset(df, contact_rate_mean==cr & `Vaccine efficacy`==ve)
        
        g <- ggplot(df, aes(x=ageGroup, y=mn, fill=factor(interv_start))) +
            geom_point(alpha=0.75, size=3, shape=22) +
            geom_point(size=3, shape=0, alpha=0.5,color='black') +
            facet_wrap(interv_cvg_rate ~  outcome, scales = 'free',ncol = 3) +
            guides(fill=guide_legend(title="Vacc. Lag (days)")) +
            #coord_cartesian(ylim=c(0,1)) + 
            ggtitle(paste('CR =',cr,'\n VaxEff =',ve))+
            scale_fill_brewer(palette = 'RdYlGn')+
            xlab("Age group") + ylab("Mean relative reduction") +
            theme_gray()
        plot(g)
    }
    
    
    pdf(file = paste0('../results/Fig_OutcomeAge_CRVE.pdf'), 
        width = 11,height = 10)
    for(i in 1:nrow(q)) subplot(q[i,1],q[i,2])
    for(i in 1:nrow(q)) subplot.2(q[i,1],q[i,2])
    dev.off()
}


mypalette <- c('-28'='red', '-14'='orange', 
               '0'='grey',
               '14'='royalblue1','28'='royalblue3')
mypalette2 <- c('1200'='red', '600'='orange', 
               '300'='royalblue1','100'='royalblue4')


figure.S1a <- function(zlist.ag) {
    
    df <- do.call('rbind.data.frame', zlist.ag)
    df$ageGroup <- factor(df$ageGroup, levels = c("0_5", "5_18", "18_65","65_over"))
    ar <- unique(df$interv_cvg_rate)
    
    # Explicit name for plot:
    df$VE <- paste('VE =',df$interv_efficacy)
    df$AG <- paste('Age group =',df$ageGroup)
    df$AG <- gsub('_',' to ',df$AG)
    df$AG <- gsub('to over','and over',df$AG)
    df$CR <- paste('[Ro] =',df$contact_rate_mean)
    df <- subset(df, interv_target=='priority_age_frailty')
    
    mypalette <- c('-28'='red', '-14'='orange', 
                   '0'='grey',
                   '14'='royalblue1','28'='royalblue3')
    
    g <- ggplot(df, aes(x=interv_cvg_rate, y=mn, color=factor(interv_start))) +
        geom_line(alpha=0.5, size=2) +
        geom_point(alpha=0.7, size=3) +
        geom_point(size=3, shape=1) +
        scale_x_continuous(breaks=ar) +
        facet_grid( VE + CR ~ AG) +
        guides(color=guide_legend(title="Vacc. Lag (days)")) +
        ggtitle('Figure S1 a')+
        scale_color_manual(values=mypalette) +
        theme(panel.grid.minor.x = element_blank()) +
        xlab("Vaccine administration rate (per 100,000 per day)") + ylab("Mean relative reduction") 
    
    pdf(file = paste0('../results/Fig_S1a.pdf'), width = 11, height = 10)    
    plot(g)
    dev.off()
}

figure.S1 <- function(zlist) {
    
    df <- do.call('rbind.data.frame', zlist)
    ar <- unique(df$interv_cvg_rate)
    
    # Explicit name for plot:
    df$VE <- paste('VE =',df$interv_efficacy)
    df$CR <- paste('[Ro] =',df$contact_rate_mean)
    df <- subset(df, interv_target=='priority_age_frailty')
    
    g <- ggplot(df, aes(x=interv_cvg_rate, y=mn, color=factor(interv_start))) +
        geom_line(alpha=0.5, size=2) +
        geom_point(alpha=0.7, size=3) +
        geom_point(size=3, shape=1) +
        scale_x_continuous(breaks=ar) +
        facet_grid( VE ~ CR ) +
        guides(color=guide_legend(title="Vacc. Lag (days)")) +
        ggtitle('Figure S1')+
        theme(panel.grid.minor.x = element_blank()) +
        scale_color_manual(values=mypalette) +
        xlab("Vaccine administration rate (per 100,000 per day)") + ylab("Mean relative reduction") 
    
    pdf(file = paste0('../results/Fig_S1.pdf'), width = 11, height = 10)    
    plot(g)
    dev.off()
}

figure.S2 <- function(zlist.ag) {
    
    df <- do.call('rbind.data.frame', zlist.ag)
    df$ageGroup <- factor(df$ageGroup, levels = c("0_5", "5_18", "18_65","65_over"))
    ar <- unique(df$interv_cvg_rate)
    
    df <- subset(df, interv_target=='priority_age_frailty')
    
    # Explicit name for plot:
    df$VE <- paste('VE =',df$interv_efficacy)
    df$AG <- paste('Age group =',df$ageGroup)
    df$AG <- gsub('_',' to ',df$AG)
    df$AG <- gsub('to over','and over',df$AG)
    df$CR <- paste('[Ro] =',df$contact_rate_mean)
    
    mypalette <- c('-28'='red', '-14'='orange', 
                   '0'='grey',
                   '14'='royalblue1','28'='royalblue3')
    dodge <- position_dodge(width=185)
    
    g <- ggplot(df, aes(x=interv_cvg_rate, y=mn, 
                        fill=factor(interv_start))) +
        geom_bar(position = 'dodge', stat = 'identity', alpha=0.7) +
        geom_errorbar(aes(ymin=qlo, ymax=qhi,
                          color = factor(interv_start)),
                      width=100,
                      size = 0.2,
                      position = dodge)+
        scale_x_continuous(breaks=ar) +
        facet_grid( VE + CR ~ AG) +
        guides(fill=guide_legend(title="Vacc. Lag (days)")) +
        ggtitle('Figure S2')+
        scale_fill_manual(values=mypalette) +
        scale_color_manual(values=mypalette) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank()) +
        xlab("Vaccine administration rate (per 100,000 per day)") + ylab("Mean relative reduction") 
    
    pdf(file = paste0('../results/Fig_S2.pdf'), width = 11, height = 10)    
    plot(g)
    dev.off()
}

figure.S3 <- function(zlist.ag){
    
    df <- do.call('rbind.data.frame', zlist.ag)
    # df$outcome  <- factor(df$outcome, levels = c("Symptomatic Infections", "Hospitalized", "Deaths"))
    df$ageGroup <- factor(df$ageGroup, levels = c("0_5", "5_18", "18_65","65_over"))
    
    df <- subset(df, interv_target=='priority_age_frailty')
    
    df$VE <- paste('VE =',df$interv_efficacy)
    df$CR <- paste('[Ro] =',df$contact_rate_mean)
    
    g <- ggplot(df, aes(x=ageGroup, y=mn, 
                        fill=factor(interv_start),
                        color=factor(interv_start))) +
        geom_line(aes(group=factor(interv_start))) +
        geom_point(alpha=0.75, size=3, shape=22) +
        facet_grid( VE + CR ~  interv_cvg_rate) +
        guides(fill=guide_legend(title="Vacc. Lag (days)"), color=FALSE) +
        scale_fill_manual(values=mypalette) + 
        scale_color_manual(values=mypalette) + 
        ggtitle('Figure S3')+
        xlab("Age group") + ylab("Mean relative reduction") +
        theme_bw()
    pdf(file = paste0('../results/Fig_S3.pdf'), width = 12,height = 8)
    plot(g)
    dev.off()
}


reformat <- function(DAT){
    DAT$VE <- paste('VE =' , DAT$interv_efficacy)
    
    # Add R0 information:
    ucr <- unique(DAT$contact_rate_mean)
    DAT$tmp <- 1.4
    DAT$tmp[DAT$contact_rate_mean==ucr[2]] <- 1.9
    DAT$R0 <- paste('Ro =',DAT$tmp)
    
    DAT$outcome <- factor(DAT$outcome, levels = c('Symptomatic Infections',
                                                  'Hospitalized',
                                                  'Deaths'))
    
    return(DAT)
}

figure.1 <- function(zlist) {
    
    DAT <- do.call('rbind.data.frame', zlist)
    DAT <- reformat(DAT)
    
    ar <- unique(DAT$interv_start)
    
    # Define and reorder VE:
    DAT$VE <- paste('VE =',DAT$interv_efficacy)
    uie <- sort(unique(DAT$interv_efficacy),decreasing = T)
    DAT$VE <- factor(x = DAT$VE, levels = paste('VE =',uie))
    
    # Select only one strategy for this plot:
    DAT <- subset(DAT,interv_target=='priority_age_frailty' )
    
    g <- ggplot(DAT, aes(x=interv_start, y=mn, 
                         fill=factor(interv_cvg_rate),
                         color=factor(interv_cvg_rate) )) +
        geom_bar(stat = 'identity', 
                 position='dodge', alpha=0.6, width=10) +
        geom_errorbar(aes(ymin=qlo,ymax=qhi),
                      size = 0.4,
                      position = 'dodge',
                      width = 10)+
        scale_x_continuous(breaks=ar) +
        facet_grid(VE ~ R0 ) +
        guides(fill = guide_legend(title="Vacc. admin. rate \n(per 100,000 per day)"),
               color = FALSE) +
        theme(panel.grid.minor.x = element_blank()) +
        theme(panel.grid.major.x = element_blank()) +
        theme(axis.text=element_text(size=18),
              axis.title=element_text(size=18,face="bold")) +
        theme(axis.title.y=element_text(margin=margin(0,20,0,0))) +
        theme(axis.title.x=element_text(margin=margin(20,0,0,0))) +
        theme(strip.text = element_text(size=18)) +
        theme(legend.text=element_text(size=16)) +
        theme(legend.title=element_text(size=16)) +
        scale_fill_manual(values=mypalette2) +
        scale_color_manual(values=mypalette2) +
        # scale_fill_brewer(palette = 'BrBG') +
        # scale_color_brewer(palette = 'BrBG') +
        xlab("Vacc. Lag (days)") + ylab("Mean relative reduction") 
    
    pdf(file = paste0('../results/figures/Figure_1.pdf'), 
        width = 12, height = 10)    
    plot(g)
    dev.off()
}


figure.1a <- function(zlist) {
    
    DAT <- do.call('rbind.data.frame', zlist)
    ar <- unique(DAT$interv_start)
    
    DAT <- subset(DAT, interv_efficacy==0.8)
    
    # Explicit name for plot:
    DAT$CR <- paste('[Ro] =',DAT$contact_rate_mean)
    DAT$it <- NA
    DAT$it[DAT$interv_target=='never_sympt'] <- 'Random'
    DAT$it[DAT$interv_target=='priority_age_frailty'] <- 'Priority'
    
    g <- ggplot(DAT, aes(x=interv_start, y=mn, 
                        fill=factor(interv_cvg_rate),
                        color=factor(interv_cvg_rate) )) +
        geom_bar(stat = 'identity', position='dodge', alpha=0.6) +
        geom_errorbar(aes(ymin=qlo,ymax=qhi),
                      position = 'dodge')+
        scale_x_continuous(breaks=ar) +
        facet_grid(it ~ CR ) +
        guides(fill = guide_legend(title="Vacc. admin. rate \n(per 100,000 per day)"),
               color = FALSE) +
        ggtitle('Figure 1')+
        theme(panel.grid.minor.x = element_blank()) +
        theme(panel.grid.major.x = element_blank()) +
        scale_fill_manual(values=mypalette2) +
        scale_color_manual(values=mypalette2) +
        # scale_fill_brewer(palette = 'BrBG') +
        # scale_color_brewer(palette = 'BrBG') +
        xlab("Vacc. Lag (days)") + ylab("Mean relative reduction") 
    
    pdf(file = paste0('../results/Fig_1a.pdf'), width = 11, height = 10)    
    plot(g)
    dev.off()
}



figure.1.b <- function(zlist) {
    
    DAT <- do.call('rbind.data.frame', zlist)
    ar <- unique(DAT$interv_start)
    DAT <- subset(DAT, interv_efficacy==0.8)
    
    # Explicit name for plot:
    DAT$CR <- paste('[Ro] =',DAT$contact_rate_mean)
    DAT$it <- NA
    DAT$it[DAT$interv_target=='never_sympt'] <- 'Random'
    DAT$it[DAT$interv_target=='priority_age_frailty'] <- 'Priority'
    
    g <- ggplot(DAT, aes(x=interv_start, y=mn, color=factor(interv_cvg_rate))) +
        geom_line(alpha=0.5, size=2) +
        geom_point(alpha=0.8, size=3) +
        # geom_point(size=3, shape=1) +
        scale_x_continuous(breaks=ar) +
        facet_grid(it ~ CR ) +
        guides(color=guide_legend(title="Vacc. admin. rate \n(per 100,000 per day)")) +
        ggtitle('Figure 1b')+
        theme(panel.grid.minor.x = element_blank()) +
        scale_color_manual(values=mypalette2) +
        # scale_color_brewer(palette = 'BrBG') +
        xlab("Vacc. Lag (days)") + ylab("Mean relative reduction") 
    
    pdf(file = paste0('../results/Fig_1b.pdf'), width = 11, height = 10)    
    plot(g)
    dev.off()
}



figure.1.c.alloutcome <- function(zlist.full) {
    
    DAT <- do.call('rbind.data.frame', zlist.full)
    ar <- unique(DAT$interv_start)
    #DAT <- subset(DAT, interv_efficacy==0.9)
    
    # Explicit name for plot:
    DAT$CR <- paste('[Ro] =',DAT$contact_rate_mean)
    DAT$it <- NA
    DAT$it[DAT$interv_target=='never_sympt'] <- 'Random'
    DAT$it[DAT$interv_target=='priority_age_frailty'] <- 'Priority'
    
    g <- ggplot(DAT, aes(x=interv_start, y=mn, color=factor(interv_cvg_rate))) +
        geom_line(alpha=0.5, size=2, aes(linetype=factor(it))) +
        geom_point(alpha=0.8, size=3) +
        scale_x_continuous(breaks=ar) +
        facet_grid( outcome ~ CR + interv_efficacy ) +
        guides(color=guide_legend(title="Vacc. admin. rate \n(per 100,000 per day)")) +
        theme(panel.grid.minor.x = element_blank()) +
        scale_color_manual(values=mypalette2) +
        xlab("Vacc. Lag (days)") + ylab("Mean relative reduction") 
    
    pdf(file = paste0('../results/Fig_1c_alloutcomes.pdf'), width = 16, height = 10)    
    plot(g)
    dev.off()
}




figure.1.c <- function(zlist) {
    
    DAT <- do.call('rbind.data.frame', zlist)
    ar <- unique(DAT$interv_start)
    DAT <- subset(DAT, interv_efficacy==0.8)
    
    # Explicit name for plot:
    DAT$CR <- paste('[Ro] =',DAT$contact_rate_mean)
    DAT$it <- NA
    DAT$it[DAT$interv_target=='never_sympt'] <- 'Random'
    DAT$it[DAT$interv_target=='priority_age_frailty'] <- 'Priority'
    
    g <- ggplot(DAT, aes(x=interv_start, y=mn, color=factor(interv_cvg_rate))) +
        geom_line(alpha=0.5, size=2, aes(linetype=factor(it))) +
        geom_point(alpha=0.8, size=3) +
        # geom_point(size=3, shape=1) +
        scale_x_continuous(breaks=ar) +
        facet_grid(~ CR ) +
        guides(color=guide_legend(title="Vacc. admin. rate \n(per 100,000 per day)")) +
        # ggtitle('Figure 1c')+
        theme(panel.grid.minor.x = element_blank()) +
        scale_color_manual(values=mypalette2) +
        # scale_color_brewer(palette = 'BrBG') +
        xlab("Vacc. Lag (days)") + ylab("Mean relative reduction") 
    
    pdf(file = paste0('../results/Fig_1c.pdf'), width = 11, height = 10)    
    plot(g)
    dev.off()
}

process.outputs <- function(df,dir,file.scen.prm.list) {
    # Retrieve scenarios definitions:
    spl <- read.csv(file.scen.prm.list) # DEBUG::  file.scen.prm.list <- 'scenario-prm-list.csv'
    names(df)[names(df)=='scenario'] <- 'scenario_id' 
    
    # Merge results and scenarios definition:
    x <- join(df, spl, by='scenario_id')
    
    z.sympt.ag <- output.stats(x, resp.var = 'rel.d.sympt',do.ageGroup = TRUE)
    z.sympt    <- output.stats(x, resp.var = 'rel.d.sympt',do.ageGroup = FALSE)
    z.hosp     <- output.stats(x, resp.var = 'rel.d.hosp', do.ageGroup = FALSE)
    z.death    <- output.stats(x, resp.var = 'rel.d.death',do.ageGroup = FALSE)
    
    titlelist <- list('Symptomatic Infections','Hospitalized','Deaths')
    
    zlist <- list(z.sympt)  #, z.hosp, z.death)
    zlist <- lapply(zlist, FUN = format.df.for.figure)
    
    zlist.ag <- list(z.sympt.ag)  #, z.hosp, z.death)
    zlist.ag <- lapply(zlist.ag, FUN = format.df.for.figure)
    
    zlist.full <- list(z.sympt, z.hosp, z.death)
    zlist.full <- lapply(zlist.full, FUN = format.df.for.figure)
    
    for(i in seq_along(zlist)) zlist[[i]]$outcome <- titlelist[[i]]
    for(i in seq_along(zlist.ag)) zlist.ag[[i]]$outcome <- titlelist[[i]]
    
    for(i in seq_along(zlist.full)) zlist.full[[i]]$outcome <- titlelist[[i]]
    
    return(list(zlist = zlist,
                zlist.ag = zlist.ag,
                zlist.full = zlist.full))
}

figures.maintext <-function(df,dir,file.scen.prm.list){
    
    X <- process.outputs(df,dir,file.scen.prm.list)
    
    zlist      <- X[['zlist']]
    zlist.ag   <- X[['zlist.ag']]
    zlist.full <- X[['zlist.full']]
    
    # Plot and save :
    figure.1(zlist)
    
    
    # figure.1.b(zlist)
    # figure.1.c(zlist)
    # 
    # figure.1.c.alloutcome(zlist.full)
    # 
    # figure.S1a(zlist.ag) 
    # figure.S1(zlist) 
    # figure.S2(zlist.ag) 
    # figure.S3(zlist.ag)
}



