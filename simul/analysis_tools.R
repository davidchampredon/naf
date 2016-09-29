
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




plot.age.contact.matrix <- function(res) {
	
	# Desired contact assortativity:
	y <- res$contactAssort
	D <- matrix(unlist(res$contactAssort), ncol = length(y))
	
	# Effective contacts from simulation:
	x <- res$wiw_ages
	m <- matrix(unlist(x), ncol=length(x))
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
	A.plot <- log(1+A)
	
	# -- PLOTS --
	
	par(mfrow = c(1,2))
	
	na <- ncol(A)
	D.plot <- D[1:(na+1),1:(na+1)]
	
	image(x=0:na, y=0:na, z=D.plot, 
		  col = topo.colors(12), 
		  main='Input contact assortativity',
		  xlab = 'age', ylab='age', las=1)
	abline(a=0,b=1,lty=2); grid()
	
	image(A.plot,x = 1:m.max, y=1:m.max, zlim = c(0,max(A.plot)), 
		  ylab = 'infector\'s age',
		  xlab = 'infectee\'s age',
		  las = 1,
		  main = 'Simulated effective contact age matrix',
		  col  = topo.colors(12))
	abline(a=0,b=1,lty=2); grid()
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

plot.sched <- function(pop){
	
	df <- ddply(pop,c('sched_type'), summarize, 
				m = mean(n_secondary_cases),
				sd = sd(n_secondary_cases))
	
	
	g <- ggplot(df) + geom_point(aes(x=factor(sched_type), y=m), size=6)
	g <- g + geom_segment(aes(x=sched_type,xend=sched_type, y=m-sd, yend=m+sd), size=2)
	g <- g + ggtitle('Mean & sd of number secondary cases')
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
	
	g.R.symptom <- ggplot(pop.transm) +geom_histogram(aes(x=n_secondary_cases, ..density..,
														  fill = factor(was_symptomatic),
														  colour = factor(was_symptomatic)),
													  position=position_dodge(width=0.8),
													  binwidth = 1,
													  alpha = 0.7)
	
	g.R.symptom <- g.R.symptom + scale_x_continuous(breaks =  seq(0,nsmax+1,by=1))
	
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
	
	g.R.treat <- ggplot(pop.transm) +geom_histogram(aes(x=n_secondary_cases, ..density..,
														fill = factor(is_treated),
														colour = factor(is_treated)),
													position=position_dodge(width=0.8),
													binwidth = 1,
													alpha = 0.7)
	
	g.R.treat <- g.R.treat + scale_x_continuous(breaks =  seq(0,nsmax+1,by=1))
	
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
	
	
	# === Schedules ===
	
	g.sched.R <- plot.sched(pop)
	
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
				 g.sympt.imm.cell.dist,
				 g.sched.R
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
	g.SR <- g.SR + geom_step(aes(y=n_vaccinated),colour='springgreen2', linetype = 1)
	g.SR <- g.SR + geom_step(aes(y=n_treated),colour='tomato', linetype = 1)
	g.SR <- g.SR + ggtitle("Susceptible, recovered, vaccinated, treated") + ylab("")
	
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




plot.sp.one <- function(pop, world.prm, name.in, name.out){
	
	undefined.id <- 999999999
	
	for(id.au in unique(pop$id_au))
	{
		pop.au <- subset(pop, id_au==id.au)
		pop.au <- pop.au[pop.au[,name.out] < undefined.id, ]
		df <- ddply(pop.au, name.out, function(x){c(sz=length(x[,name.out]))})
		df <- ddply(df,'sz',summarise, n=length(sz))
		
		sz.out <- df$sz
		pr.out <- df$n/sum(df$n)
		sz.in <- world.prm[[name.in]][[id.au+1]]
		pr.in <- world.prm[[paste0(name.in,'_proba')]][[id.au+1]]
		
		plot(x = sz.out, 
			 y = pr.out, 
			 main = paste0('Size Distribution in AU_',id.au),
			 ylim = range(pr.out, pr.in),
			 xlim = range(sz.in, sz.out),
			 xlab = name.in, ylab='prop', las=1,
			 typ ='l', lwd = 2)
		
		lines(x = sz.in, 
			   y = pr.in, 
			   pch = 3, cex=2, lwd=2,
			   col='red', typ='o')
		points(x = sz.out, y = pr.out, 
			   pch=1, lwd = 2, cex=1.5)
		
		legend(x='topright', pch=c(3,1), col=c('red','black'), 
			   pt.lwd = 3,pt.cex=1,
			   legend=c('Target','Realized'))
	}
}



plot.sp.sz.distrib <- function(pop,world.prm) {
	
	name.in.vec  <- c('hh_size','wrk_size', 'school_size', 'pubt_size', 'other_size')
	name.out.vec <- c('id_hh',  'id_wrk',   'id_school',   'id_pubTr',  'id_other')
	
	par(mfrow=c(5,2), cex.lab=2, cex.axis=2)
	for(i in seq_along(name.out.vec)){
		plot.sp.one(pop, world.prm, 
					name.in  = name.in.vec[i], 
					name.out = name.out.vec[i])
	}
}


plot.share.same.hh <- function(pop) {
	# Proportion of indiv from the _same_ household
	# sharing the same social places (other than household!)
	
	same.hh.one <- function(sptype) {
		pop$key <- paste(pop$id_hh, pop[,sptype], sep='_')
		x <-ddply(pop, sptype, function(x){ c(  q = length(unique(x[,'key'])), 
												n = length(x[,sptype])   )})
		x <- subset(x, x[,sptype]<999999999)
		x$ratio <- 1-x$q/x$n
		return(data.frame(sptype,ratio=x$ratio))
	}
	sptype.vec <- c('id_wrk','id_pubTr','id_other','id_school')
	y <- lapply(sptype.vec, same.hh.one)
	yy <- do.call('rbind',y)
	par(mfrow=c(1,1), cex.lab=2, cex.axis=2)
	boxplot(ratio~sptype, data=yy, ylim=c(0,1), col='lightgrey',
			main='Proportion of indiv sharing the same household\n(all AU combined)')
	grid()
}


