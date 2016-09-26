
library(plyr)
library(gridExtra)


synthetic_age_adult <- function(age.adult){
	# Create a synthetic age distribution for adults.

	age.thres <- 60
	age.max <- max(age.adult)
	idx <- which(age.adult > age.thres)
	# Decline of older adults:
	rel.prop <- (age.max - age.adult[idx])/(age.max - age.thres)
	
	p.adult <- c( rep(1.0, idx[1]-1) , rel.prop)
	p.adult <- p.adult/sum(p.adult)
	return(p.adult)
}


plot.binomial.regression <- function(dat, xvar, binomial_response, title) {
	n <- nrow(dat)
	
	g <- ggplot(dat) 
	if(n<1000) g <- g + geom_point(aes_string(x=xvar,y=binomial_response), alpha=0.3) 
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


plot.age.contact.matrix <- function(x) {
	# x <- res$wiw_ages
	m <- matrix(unlist(x), ncol=length(x))
	
	hist(m[,1],breaks=100, ylim=c(0,2000))
	
	m.max <- ceiling(max(m))
	
	A <- matrix(nrow = m.max, ncol=m.max, data = 0)
	ages <- 1:m.max
	row.names(A) <- ages
	colnames(A)  <- ages
	
	for(q in 1:nrow(m)){
		i <- ceiling(m[q,1])
		j <- ceiling(m[q,2])
		A[i,j] <- A[i,j] + 1
	}
	A <- sqrt(A) #log(A+1)
	
	image(A,x = 1:m.max, y=1:m.max, zlim = c(0,max(A)), 
		  ylab = 'infector\'s age',
		  xlab = 'infectee\'s age',
		  main = 'Age-contact matrix',
		  col  = topo.colors(12))
}


plot.n.contacts <- function(nc){
	df0 <- data.frame(time=nc$time, uid=nc$uid, n=nc$nContacts)
	df0$timeround <- ceiling(df0$time)
	
	df    <- ddply(df0, c('timeround','uid'),summarize, ncontacts = sum(n))
	df.ts <- ddply(df0, c('timeround'),summarize, tot.contacts = sum(n))
	
	gts <- ggplot(df.ts, aes(x=timeround, y=tot.contacts))+ geom_step()
	gts <- gts + scale_y_log10() + ggtitle('Total number of contacts')+xlab('time')+ylab('')
	
	# Distribution
	m <- mean(df$ncontacts)
	ci <- 0.95
	qlo <- quantile(df$ncontacts, probs=(1-ci)/2)
	qhi <- quantile(df$ncontacts, probs=0.5 + ci/2)
	
	g <- ggplot(df,aes(x=ncontacts)) + geom_histogram(binwidth=1,fill='gold2',colour='gold3') 
	g <- g + ggtitle(paste0('Distribution for the number of contacts (mean = ',round(m,2),' ; ',
							ci*100,'%CI: ',round(qlo,2),' -- ',round(qhi,2),')'))
	g <- g + xlab('Contacts per day, per individual')#+scale_y_log10()
	g <- g + geom_vline(xintercept=m, linetype=2, colour='gold3', size=2)
	
	grid.arrange(gts,g)	
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
	
	g.death.imm.hum <- plot.binomial.regression(dat = pop, 
											xvar='immunity_hum',
											binomial_response = 'is_alive',
											title='Death and Humoral Immunity')
	
	g.death.imm.cell <- plot.binomial.regression(dat = pop, 
											xvar='immunity_cell',
											binomial_response = 'is_alive',
											title='Death and Cellular Immunity')
	
	
	g.death.frailty.dist <- plot.density.categ(dat = pop, 
											   xvar='frailty',
											   categ = 'is_alive',
											   title='Frailty distribution (by survival)')
	
	g.death.imm.hum.dist <- plot.density.categ(dat = pop, 
											   xvar='immunity_hum',
											   categ = 'is_alive',
											   title='Humoral Immunity distribution (by survival)')
	
	g.death.imm.cell.dist <- plot.density.categ(dat = pop, 
												xvar='immunity_cell',
												categ = 'is_alive',
												title='Cellular Immunity distribution (by survival)')
	
	
	
	
	# ==== Immunity and frailty ====
	
	pop$age.round <- round(pop$age,0)
	
	age.imm.fra <- ddply(pop,c('age.round'),summarize, 
						 imm.hum = mean(immunity_hum), 
						 imm.cell = mean(immunity_cell), 
						 fra = mean(frailty))
	age.imm.fra
	
	g.age.imm.hum <- ggplot(age.imm.fra)+geom_line(aes(x=age.round,y=imm.hum)) +coord_cartesian(ylim = c(0,1))
	g.age.imm.hum <- g.age.imm.hum + ggtitle('Humoral Immunity by age') + xlab('age') + ylab('immunity')
	g.age.imm.hum <- g.age.imm.hum + geom_smooth(aes(x=age.round,y=imm.hum), colour='red3',size=2,se = F)
	
	g.age.imm.cell <- ggplot(age.imm.fra)+geom_line(aes(x=age.round,y=imm.cell)) +coord_cartesian(ylim = c(0,1))
	g.age.imm.cell <- g.age.imm.cell + ggtitle('Cellular Immunity by age') + xlab('age') + ylab('immunity')
	g.age.imm.cell <- g.age.imm.cell + geom_smooth(aes(x=age.round,y=imm.cell), colour='red3',size=2,se = F)
	
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
	
	g.imm.hum.hosp <- plot.binomial.regression(dat = pop, xvar = 'immunity_hum',
											 binomial_response = 'hosp',
											 title = 'Hospitalization and humoral immunity')
	
	g.imm.cell.hosp <- plot.binomial.regression(dat = pop, xvar = 'immunity_cell',
											   binomial_response = 'hosp',
											   title = 'Hospitalization and cellular immunity')
	
	
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

	g.sympt.imm.hum.dist <- plot.density.categ(dat = pop, 
											   xvar ='immunity_hum',
											   categ = 'was_symptomatic',
											   title ='Humoral Immunity distribution\n(by symptomatic status)')
	
	g.sympt.imm.cell.dist <- plot.density.categ(dat = pop, 
												xvar ='immunity_cell',
												categ = 'was_symptomatic',
												title ='Cellular Immunity distribution\n(by symptomatic status)')
	
	
	sympt.vax <- ddply(pop.inf,c("is_vaccinated","was_symptomatic"), summarize, n=length(id_indiv))
	sympt.vax
	
	g.sympt.vax <- ggplot(sympt.vax)+geom_bar(aes(x=factor(is_vaccinated), y=n, fill=factor(was_symptomatic)),
												  stat='identity', position='dodge')
	g.sympt.vax <- g.sympt.vax + ggtitle('Symptomatic infection and vaccination')
	

	### ==== Final ====
	
	grid.arrange(g.age, 
				 g.age.imm.hum,
				 g.age.imm.cell,
				 g.age.fra,
				 g.age.death.dist,
				 g.death.frailty,
				 g.death.imm.hum,
				 g.death.imm.cell,
				 g.death.frailty.dist,
				 g.death.imm.hum.dist,
				 g.death.imm.cell.dist,
				 g.imm.hum.hosp,
				 g.imm.cell.hosp,
				 g.frail.hosp,
				 g.age.hosp, 
				 g.treat.hosp,
				 g.vax.hosp,
				 g.sympt.vax,
				 g.sympt.imm.hum.dist,
				 g.sympt.imm.cell.dist
				)
	
	grid.arrange( g.dol.drawn, 
				  g.doi.drawn,
				  g.dobh.drawn,
				  g.doh.drawn, g.R,
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

plot.ts.sp <- function(ts, facets=FALSE){

	zz <- data.frame(matrix(unlist(ts),ncol = length(ts)))
	names(zz) <- names(ts)
	
	myconv <- function(x) {
		return(as.numeric(as.character(x)))
	}
	zz$time <- myconv(zz$time)
	zz$id_sp<- myconv(zz$id_sp)
	zz$nS <- myconv(zz$nS)
	zz$nE <- myconv(zz$nE)
	
	zz$timeround <- floor(zz$time)
	df <- ddply(zz,c('timeround','type'),summarize, n=sum(nE))
	
	g <- ggplot(df,aes(x=timeround, y=n, colour=type, shape=type))
	g <- g + geom_point()+geom_line(size=0.2)
	if(facets) g <- g +facet_wrap(~type)
	g <- g + ggtitle("Exposed individuals, by SP types") + xlab('day')
	g <- g + scale_y_log10()
	return(g)
}
