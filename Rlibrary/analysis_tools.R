
library(plyr)
library(gridExtra)


synthetic_age_adult <- function(age.adult){
	# Create a synthetic age distribution for adults.

	age.thres <- 60
	age.max <- 90
	idx <- which(age.adult > age.thres)
	# Decline of older adults:
	rel.prop <- (age.max - age.adult[idx])/(age.max - age.thres)
	
	p.adult <- c( rep(1.0, idx[1]-1) , rel.prop)
	p.adult <- p.adult/sum(p.adult)
	return(p.adult)
}


plot.binomial.regression <- function(dat, xvar, binomial_response, title) {
	g <- ggplot(dat) + geom_point(aes_string(x=xvar,y=binomial_response), alpha=0.3) 
	g <- g + geom_smooth(aes_string(x=xvar, y=binomial_response), 
						 method = "glm", 
						 method.args = list(family = "binomial"), 
						 colour='red3',size=2,se = F)
	g <- g + ggtitle(title)
	return(g)
}

plot.density.categ <- function(dat, xvar, categ, title) {
	
	dat[,categ] <- as.factor(dat[,categ])
	g <- ggplot(dat) 
	
	if (categ != ''){
		g <- g + geom_density(aes_string(x=xvar,
										 fill=categ,
										 colour=categ),
							  alpha=0.3)
	}
	g <- g + geom_line(stat='density',
					   aes_string(x=xvar),
					   size = 1.5)
	g <- g + ggtitle(title)
	return(g)
}


plot.population <- function(pop) {
	
	pop$hosp <- as.numeric( as.logical(pop$is_discharged+pop$is_hosp) )
	pop.hosp <- subset(pop, hosp>0)
	pop.transm <- subset(pop, n_secondary_cases>0)
	pop.inf <- subset(pop, doi_drawn>0)
	
	# ==== Age ====
	
	g <- ggplot(pop)
	g.age <- g + geom_histogram(aes(x=age), binwidth=2, fill='darkgrey', colour = 'black')
	g.age <- g.age + ggtitle('Age distribution')
	
	g.age.death.dist <- plot.density.categ(dat = pop, 
										   xvar='age',
										   categ = 'is_alive',
										   title='Age distribution (by survival)')
	
	# ==== Death ====
	
	
	g.death.frailty <- plot.binomial.regression(dat = pop, 
												xvar='frailty',
												binomial_response = 'is_alive',
												title='Death and Frailty')
	
	g.death.imm <- plot.binomial.regression(dat = pop, 
											xvar='immunity',
											binomial_response = 'is_alive',
											title='Death and Immunity')
	
	g.death.frailty.dist <- plot.density.categ(dat = pop, 
											   xvar='frailty',
											   categ = 'is_alive',
											   title='Frailty distribution (by survival)')
	g.death.imm.dist <- plot.density.categ(dat = pop, 
										   xvar='immunity',
										   categ = 'is_alive',
										   title='Immunity distribution (by survival)')
	
	
	
	# ==== Immunity and frailty ====
	
	pop$age.round <- round(pop$age,0)
	
	age.imm.fra <- ddply(pop,c('age.round'),summarize, imm = mean(immunity), fra = mean(frailty))
	age.imm.fra
	
	g.age.imm <- ggplot(age.imm.fra)+geom_line(aes(x=age.round,y=imm)) +coord_cartesian(ylim = c(0,1))
	g.age.imm <- g.age.imm + ggtitle('Immunity by age') + xlab('age') + ylab('immunity')
	g.age.imm <- g.age.imm + geom_smooth(aes(x=age.round,y=imm), colour='red3',size=2,se = F)
	
	g.age.fra <- ggplot(age.imm.fra)+geom_line(aes(x=age.round,y=fra)) +coord_cartesian(ylim = c(0,1))
	g.age.fra <- g.age.fra + ggtitle('Frailty by age') + xlab('age') + ylab('frailty')
	g.age.fra <- g.age.fra + geom_smooth(aes(x=age.round,y=fra), colour='red3',size=2,se = F)
	
	
	# ==== DOI & DOL ====
	
	linew = 1
	alpha = 0.8
	m_dol_drawn = mean(pop.inf$dol_drawn)
	m_doi_drawn = mean(pop.inf$doi_drawn)

	g <- ggplot(pop.inf)
	
	g.dol.drawn <- g + geom_histogram(aes(x=dol_drawn), fill='orange',colour='orange', size=linew, alpha=alpha, binwidth = 0.5)
	g.dol.drawn <- g.dol.drawn + geom_vline(xintercept = m_dol_drawn, colour='orange',linetype = 2)
	g.dol.drawn <- g.dol.drawn + ggtitle(paste0("DOL drawn distribution, among infected (mean = ", round(m_dol_drawn,2),")"))
	
	g.doi.drawn <- g + geom_histogram(aes(x=doi_drawn), fill='red',colour='red', size=linew, alpha=alpha, binwidth = 0.5)
	g.doi.drawn <- g.doi.drawn + geom_vline(xintercept = m_doi_drawn, colour='red',linetype = 2)
	g.doi.drawn <- g.doi.drawn + ggtitle(paste0("DOI drawn distribution, among infected (mean = ", round(m_doi_drawn,2),")"))
	
	
	m_dobh_drawn = mean(pop.hosp$dobh_drawn)
	m_doh_drawn  = mean(pop.hosp$doh_drawn)
	
	g.hosp <- ggplot(pop.hosp)
	g.dobh.drawn <- g.hosp + geom_histogram(aes(x=dobh_drawn), fill='purple',colour='purple', size=linew, alpha=alpha, binwidth = 0.5)
	g.dobh.drawn <- g.dobh.drawn + geom_vline(xintercept = m_dobh_drawn, colour='purple',linetype = 2)
	g.dobh.drawn <- g.dobh.drawn + ggtitle(paste0("DOBH drawn distribution, among hospitalized (mean = ", round(m_dobh_drawn,3),")"))
	
	g.doh.drawn <- g.hosp + geom_histogram(aes(x=doh_drawn), fill='green4',colour='green4', size=linew, alpha=alpha, binwidth = 0.5)
	g.doh.drawn <- g.doh.drawn + geom_vline(xintercept = m_doh_drawn, colour='green4',linetype = 2)
	g.doh.drawn <- g.doh.drawn + ggtitle(paste0("DOH drawn distribution, among hospitalized (mean = ", round(m_doh_drawn,3),")"))
	
	
	### ==== Hospitalizations ==== 
	
	g.frail.hosp <- plot.binomial.regression(dat = pop, xvar = 'frailty',
											 binomial_response = 'hosp',
											 title = 'Hospitalization and frailty')
	
	g.imm.hosp <- plot.binomial.regression(dat = pop, xvar = 'immunity',
											 binomial_response = 'hosp',
											 title = 'Hospitalization and immunity')
	
	g.age.hosp <- plot.binomial.regression(dat = pop, xvar = 'age',
										   binomial_response = 'hosp',
										   title = 'Hospitalization and age')

	
	pop.treat.hosp <- ddply(pop,c('hosp','is_treated'),summarize, n=length(id_indiv))
	g.treat.hosp <- ggplot(pop.treat.hosp) + geom_bar(aes(x = factor(hosp), 
														  y = n, 
														  fill = factor(is_treated)),
													  position = 'dodge', stat = 'identity')
	g.treat.hosp <- g.treat.hosp + ggtitle('Hospitalization and treatment')
	
	
	pop.vax.hosp <- ddply(pop,c('hosp','is_vaccinated'),summarize, n=length(id_indiv))
	g.vax.hosp <- ggplot(pop.vax.hosp) + geom_bar(aes(x = factor(hosp), 
													  y = n, 
													  fill = factor(is_vaccinated)),
												  position = 'dodge', stat = 'identity')
	g.vax.hosp <- g.vax.hosp + ggtitle('Hospitalization and vaccination')
	
	
	### ==== Secondary cases distribution ==== 
	
	R0 <- mean(pop.transm$n_secondary_cases)
	nsmax <- max(pop.transm$n_secondary_cases)
	g.R <- ggplot(pop.transm) + geom_histogram(aes(n_secondary_cases), 
											   fill = 'red3',
											   colour = 'red4',
											   breaks = seq(0,nsmax+1,by=1))
	g.R <- g.R + scale_x_continuous(breaks =  seq(0,nsmax+1,by=1))
	g.R <- g.R + geom_vline(aes(xintercept=R0), linetype = 2, size = 2)
	g.R <- g.R + ggtitle(paste0("Secondary cases distribution (R0=",
								round(R0,2),")"))
	

	# R0 symptomatic or not
	pop.transm.sum <- ddply(pop.transm,"was_symptomatic",
							summarize,
							m = mean(n_secondary_cases),
							q.lo = quantile(n_secondary_cases,probs = 0.5-0.80/2),
							q.hi = quantile(n_secondary_cases,probs = 0.5+0.80/2))
	
	pop.transm.sum
	
	g.R.symptom <- ggplot(pop.transm) +geom_density(aes(x=n_secondary_cases,
														fill = factor(was_symptomatic),
														colour = factor(was_symptomatic)),
													alpha = 0.2)
	g.R.symptom <- g.R.symptom + geom_vline(data = pop.transm.sum, 
										aes(xintercept=m, colour=factor(was_symptomatic)),
										linetype=2)
	g.R.symptom <- g.R.symptom + ggtitle(paste0("Secondary cases distribution\n (R0_asympt=",
												round(pop.transm.sum$m[1],2),
												" ; R0_sympt=",
												round(pop.transm.sum$m[2],2),
												")"))
	
	# R0 treated or not
	
	pop.transm.sum2 <- ddply(pop.transm,"is_treated",
							summarize,
							m = mean(n_secondary_cases),
							q.lo = quantile(n_secondary_cases,probs = 0.5-0.80/2),
							q.hi = quantile(n_secondary_cases,probs = 0.5+0.80/2))
	
	pop.transm.sum2
	
	g.R.treat <- ggplot(pop.transm) +geom_density(aes(x=n_secondary_cases,
													   fill = factor(is_treated),
													   colour = factor(is_treated)),
												   alpha = 0.2)
	g.R.treat <- g.R.treat + geom_vline(data = pop.transm.sum2, 
										aes(xintercept=m, colour=factor(is_treated)),
										linetype=2)
	g.R.treat <- g.R.treat + ggtitle(paste0("Secondary cases distribution\n (R0_untreated=",
											round(pop.transm.sum2$m[1],2),
											" ; R0_treated=",
											round(pop.transm.sum2$m[2],2),
											")"))

	
	### ==== Generation interval ====
	
	pop2 <- subset(pop,gi_bck>0)
	mean.gibck <- mean(pop2$gi_bck)
	sd.gibck <- sd(pop2$gi_bck)
	g.gibck <- ggplot(pop2) + geom_histogram(aes(gi_bck),
											   fill = 'skyblue1',
											   colour = 'skyblue3',
											   binwidth = 1)
	g.gibck <- g.gibck + geom_vline(xintercept=mean.gibck, 
									linetype = 2, size = 2)
	g.gibck <- g.gibck + ggtitle(paste0("Backward GI (mean=",
										round(mean.gibck,2),
										" ; sd=",
										round(sd.gibck),")"))
	
	pop2.sum <- ddply(pop2,'was_symptomatic',summarize,
					  m = mean(gi_bck))
	
	g.gibck.sympt <- ggplot(pop2) + geom_density(aes(gi_bck,
												 fill = factor(was_symptomatic),
												 colour = factor(was_symptomatic)),
												 alpha = 0.2)
	g.gibck.sympt <- g.gibck.sympt + geom_vline(data = pop2.sum , 
								aes(xintercept=m,
									colour = factor(was_symptomatic)), 
								linetype=2)
	
	g.gibck.sympt <- g.gibck.sympt + ggtitle(paste0("Backward GI (mean_asympt=",
													round(pop2.sum$m[1],2),
													" ; mean_sympt=",
													round(pop2.sum$m[2],2),
													")"))
	
	### ==== Symptomatics ====
	
	sympt.vax <- ddply(pop.inf,c("is_vaccinated","was_symptomatic"), summarize, n=length(id_indiv))
	sympt.vax
	
	g.sympt.vax <- ggplot(sympt.vax)+geom_bar(aes(x=factor(is_vaccinated), y=n, fill=factor(was_symptomatic)),
												  stat='identity', position='dodge')
	g.sympt.vax <- g.sympt.vax + ggtitle('Symptomatic infection and vaccination')
	
	
	### ==== Final ====
	
	grid.arrange(g.age, 
				 g.age.imm,
				 g.age.fra,
				 g.age.death.dist,
				 g.death.frailty,
				 g.death.imm,
				 g.death.frailty.dist,
				 g.death.imm.dist,
				 g.dol.drawn, 
				 g.doi.drawn,
				 g.dobh.drawn,
				 g.doh.drawn,
				 g.imm.hosp,
				 g.frail.hosp,
				 g.age.hosp, 
				 g.treat.hosp,
				 g.vax.hosp,
				 g.sympt.vax,
				 g.R,
				 g.R.symptom,
				 g.R.treat,
				 g.gibck,
				 g.gibck.sympt)
	
}

plot.epi.timeseries <- function(ts){
	### Plot time series
	
	ts$death_incidence <- c(ts$nD[1],diff(ts$nD))
	
	g.SR <- ggplot(ts, aes(x=time))
	g.SR <- g.SR + geom_step(aes(y=nS),colour='springgreen3', size=2) 
	g.SR <- g.SR + geom_step(aes(y=nR),colour='blue', size=2)
	g.SR <- g.SR + geom_step(aes(y=n_vaccinated),colour='springgreen2', linetype = 2)
	g.SR <- g.SR + ggtitle("Susceptible, recovered, vaccinated") + ylab("")
	
	g.inf <- ggplot(ts, aes(x=time))
	g.inf <- g.inf + geom_step(aes(y=nE),colour='orange',size=2) 
	g.inf <- g.inf + geom_step(aes(y=nIs),colour='red',size=2)
	g.inf <- g.inf + geom_step(aes(y=nIa),colour='red')
	g.inf <- g.inf + geom_step(aes(y=nIa+nIs),colour='red',linetype=3)
	g.inf <- g.inf + geom_step(aes(y=nE+nIa+nIs),colour='grey70',linetype=3)
	g.inf <- g.inf + geom_step(aes(y=n_treated), colour='cyan', linetype = 2)
	g.inf <- g.inf + geom_step(aes(y=n_vaccinated), colour='springgreen2', linetype = 2)
	g.inf <- g.inf + geom_step(aes(y=death_incidence), colour='black', linetype = 1)
	g.inf <- g.inf + ggtitle("All latent (nE, orange), Infectious (Ia & Is[bold])\n treated, vaccinated (dashed) and death (black)") + ylab("")
	
	# incidence 
	
	g.inc <- ggplot(ts, aes(x=time))
	g.inc <- g.inc + geom_step(aes(y=incidence), colour='grey') + geom_point(aes(y=incidence),size=1)
	g.inc <- g.inc + ggtitle("Raw Incidence") + ylab("")
	
	ts$cuminc <- cumsum(ts$incidence)
	g.cuminc <- ggplot(ts, aes(x=time))
	g.cuminc <- g.cuminc + geom_step(aes(y=cuminc)) 
	g.cuminc <- g.cuminc + ggtitle("Cumulative Incidence") + ylab("")
	
	ts$timeround <- floor(ts$time)
	ts2 <- ddply(ts, c('timeround'), summarize,
				 dailyinc=sum(incidence))
	g.dailyinc <- ggplot(ts2, aes(x=timeround))
	g.dailyinc <- g.dailyinc + geom_step(aes(y=dailyinc))
	# g.dailyinc <- g.dailyinc + geom_point(aes(y=dailyinc))
	g.dailyinc <- g.dailyinc + ggtitle('Daily incidence') + xlab('day')+ylab('')

	g.dailyinc.log <- g.dailyinc + scale_y_log10()
	
	g.prev <- ggplot(ts, aes(x=time))
	g.prev <- g.prev + geom_step(aes(y=prevalence))
	g.prev <- g.prev + ggtitle("Prevalence") + ylab("")
	g.prev.log <- g.prev + scale_y_log10()
	
	g.death <- ggplot(ts, aes(x=time))
	g.death <- g.death + geom_step(aes(y=nD))
	g.death <- g.death + ggtitle("Cumulative Deaths") + ylab("")
	
	g.death.inc <- ggplot(ts, aes(x=time))
	g.death.inc <- g.death.inc + geom_step(aes(y=death_incidence))
	g.death.inc <- g.death.inc + ggtitle("Deaths incidence") + ylab("")
	
	
	grid.arrange(g.SR ,
				 g.inf, 
				 g.inc,
				 g.cuminc,
				 g.dailyinc,
				 g.dailyinc.log,
				 g.prev,
				 g.prev.log,
				 g.death,
				 g.death.inc)
}
