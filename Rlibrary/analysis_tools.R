
library(plyr)
library(gridExtra)


plot.population <- function(pop) {
	
	pop$hosp <- as.numeric( as.logical(pop$is_discharged+pop$is_hosp) )
	
	g <- ggplot(pop)
	g.age <- g + geom_histogram(aes(x=age), binwidth=2, fill='darkgrey', colour = 'black')
	g.age <- g.age + ggtitle('Age distribution')
	
	linew = 1
	alpha = 0.8
	m_dol_drawn = mean(pop$dol_drawn)
	m_doi_drawn = mean(pop$doi_drawn)

	
	g.dol.drawn <- g + geom_histogram(aes(x=dol_drawn), fill='orange',colour='orange', size=linew, alpha=alpha, binwidth = 0.5)
	g.dol.drawn <- g.dol.drawn + geom_vline(xintercept = m_dol_drawn, colour='orange',linetype = 2)
	g.dol.drawn <- g.dol.drawn + ggtitle(paste0("DOL drawn distribution (mean = ", round(m_dol_drawn,2),")"))
	
	g.doi.drawn <- g + geom_histogram(aes(x=doi_drawn), fill='red',colour='red', size=linew, alpha=alpha, binwidth = 0.5)
	g.doi.drawn <- g.doi.drawn + geom_vline(xintercept = m_doi_drawn, colour='red',linetype = 2)
	g.doi.drawn <- g.doi.drawn + ggtitle(paste0("DOI drawn distribution (mean = ", round(m_doi_drawn,2),")"))
	
	pop.hosp <- subset(pop, hosp>0)
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
	
	g.frail.hosp <- ggplot(pop)
	g.frail.hosp <- g.frail.hosp + geom_smooth(aes(x=frailty, y=hosp),method = "glm", 
						 method.args = list(family = "binomial"), 
						 se = FALSE, 
						 colour='brown', size=2)
	g.frail.hosp <- g.frail.hosp + geom_point(aes(x=frailty, y=hosp), size=0.5)
	g.frail.hosp <- g.frail.hosp + ggtitle('Hospitalization and frailty')
	
	
	g.imm.hosp <- ggplot(pop)
	g.imm.hosp <- g.imm.hosp + geom_smooth(aes(x=immunity, y=hosp),method = "glm", 
											   method.args = list(family = "binomial"), 
											   se = FALSE, 
											   colour='red', size=2)
	g.imm.hosp <- g.imm.hosp + geom_point(aes(x=immunity, y=hosp), size=0.5)
	g.imm.hosp <- g.imm.hosp + ggtitle('Hospitalization and immunity')
	
	g.age.hosp <- ggplot(pop)
	g.age.hosp <- g.age.hosp + geom_smooth(aes(x=age, y=hosp),method = "glm", 
										   method.args = list(family = "binomial"), 
										   se = FALSE, 
										   colour='blue', size=2)
	g.age.hosp <- g.age.hosp + geom_point(aes(x=age, y=hosp), size=0.5)
	g.age.hosp <- g.age.hosp + ggtitle('Hospitalization and age')
	
	
	### ==== Secondary cases distribution ==== 
	
	pop.transm <- subset(pop, n_secondary_cases>0)
	R0 <- mean(pop.transm$n_secondary_cases)
	
	g.R <- ggplot(pop.transm) + geom_histogram(aes(n_secondary_cases), 
											   fill = 'red3',
											   colour = 'red4',
											   binwidth = 1)
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
	
	plot(g.R.treat)
	
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
	
	### ==== Final ====
	
	grid.arrange(g.age, 
				 g.dol.drawn, 
				 g.doi.drawn,
				 g.dobh.drawn,
				 g.doh.drawn,
				 g.imm.hosp,
				 g.frail.hosp,
				 g.age.hosp, 
				 g.R,
				 g.R.symptom,
				 g.R.treat,
				 g.gibck,
				 g.gibck.sympt)
	
}

plot.epi.timeseries <- function(ts){
	### Plot time series
	
	g.SR <- ggplot(ts, aes(x=time))
	g.SR <- g.SR + geom_step(aes(y=nS),colour='springgreen3') 
	g.SR <- g.SR + geom_step(aes(y=nR),colour='blue')
	g.SR <- g.SR + ggtitle("Susceptible and recovered") + ylab("")
	
	g.inf <- ggplot(ts, aes(x=time))
	g.inf <- g.inf + geom_step(aes(y=nE),colour='orange') 
	g.inf <- g.inf + geom_step(aes(y=nIs),colour='red',size=2)
	g.inf <- g.inf + geom_step(aes(y=nIa),colour='red')
	g.inf <- g.inf + geom_step(aes(y=n_treated),colour='green3')
	g.inf <- g.inf + ggtitle("All latent (nE, orange) and Infectious (Ia & Is[bold])") + ylab("")
	
	g.prev <- ggplot(ts, aes(x=time))
	g.prev <- g.prev + geom_step(aes(y=prevalence))
	g.prev <- g.prev + ggtitle("Prevalence") + ylab("")
	
	grid.arrange(g.SR ,
				 g.inf, 
				 g.prev)
}
