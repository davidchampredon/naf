library(profvis)
library(snowfall)
library(ggplot2)
library(plyr)
library(gridExtra)
library(tidyr)
library(parallel)


###   =====   H E L P E R   F U N C T I O N S =====



filter.out.fizzle <- function(pop.all.mc){
	### FILTER OUT THE MC ITERATIONS THAT PRODUCED FIZZLES
	
	idx.mc         <- unique(pop.all.mc$mc)
	fizz.mc        <- identify.fizzle(pop.all.mc)
	fizz.mc.idx    <- which(fizz.mc==TRUE)
	
	if(length(fizz.mc.idx)==0) 
		idx.mc.no.fizz <- idx.mc
	if(length(fizz.mc.idx)>0) 
		idx.mc.no.fizz <- idx.mc[-fizz.mc.idx]
	
	pop.nofizz <- subset(pop.all.mc, mc==idx.mc.no.fizz)
	return(pop.nofizz)
}

calc.R0 <- function(x, time.init = 3) {
	
	# 'time.init' is the maximum infection time where R0 is calculated
	# from the initial infectious individuals.
	
	# Retrieve the times of first intections
	tinf <- sort(unique(x$t_infected))
	tinf.init <- tinf[tinf < time.init]
	# Crop simulations at these early times.
	# Note that simulations are not filtered out
	# from the fizzles. (we don't care what happens later!)
	x.init <- subset(x, t_infected %in% tinf.init)
	
	R0 <- mean(x.init$n_secondary_cases)
	return(R0)
}

calc.R0.SIR <- function(pop.all.mc, 
						res.list.0,
						t.max.fit,
						do.plot = FALSE) {
	
	### ESTIMATE A R0 IMPLIED IN A SIR MODEL
	### - estimate for each MC iteration
	### - R0 = mean(estimate for each MC iteration)
	###
	
	idx.mc         <- unique(pop.all.mc$mc)
	fizz.mc        <- identify.fizzle(pop.all.mc)
	fizz.mc.idx    <- which(fizz.mc==TRUE)
	if(length(fizz.mc.idx)==0)  idx.mc.no.fizz <- idx.mc
	if(length(fizz.mc.idx)>0)   idx.mc.no.fizz <- idx.mc[-fizz.mc.idx]
	if(length(fizz.mc.idx)==length(idx.mc)) return(NA)
	
	# Filter out the fizzles
	pop.nofizz <- subset(pop.all.mc, mc %in% idx.mc.no.fizz)
	mc.vec <- unique(pop.nofizz$mc)
	tmp <- subset(pop.nofizz, gi_bck>0)
	gi_bck.mean <- mean(tmp$gi_bck)
	
	tmp2 <- ddply(.data = tmp, c('mc'), summarize, x=mean(gi_bck))
	gi.bck.mean.mc <- tmp2$x
	
	# Merging time series without fizzles:
	res.no.fizz <- list()
	for(i in seq_along(idx.mc.no.fizz)) res.no.fizz[[i]] <- res.list.0[[i]]
	ts <- merge.ts.mc(res.no.fizz, n.cpu = n.cpu)
	
	# Averaged time series:
	ts.avg <- ddply(ts,c('time'),summarize, 
					inc.m = mean(incidence),
					inc.lo = quantile(incidence,probs = 0.10),
					inc.hi = quantile(incidence,probs = 0.90))
	
	# Just the initial times, for fitting purposes:
	ts.init <- subset(ts, time < t.max.fit)
	ts.avg.init <- subset(ts.avg, time < t.max.fit)
	
	estim.R0.impl.one.mc <- function(i, t.max.fit, ts, gi.bck.mean.mc) {
		
		# Just the ith MC iteration:
		ts.i <- subset(ts, mc==i)
		
		# Find the time of first case:	
		xx <- ts.i$time
		yy <- log(ts.i$inc+1)
		ss <- cumsum(yy>0)
		idx.startepi <- max(which(ss==0))
		
		# Remove time series before start:
		xx <- xx[idx.startepi:length(xx)]
		yy <- yy[idx.startepi:length(yy)]
		
		# Just initial times:
		idx.init <- xx < t.max.fit
		xx <- xx[idx.init]
		yy <- yy[idx.init]
		
		z <- lm(formula = yy ~ xx)
		# plot(x=xx, y=yy)
		# abline(a=z$coefficients[1], b=z$coefficients[2])
		
		# r estimation on the averaged time-series:
		r.i <-  z$coefficients[2]
		
		# R0 based on SIR-like formula:
		R0.i<- 1 + r.i * gi.bck.mean.mc[i]
		return(list(R0=R0.i, r=r.i, i0=z$coefficients[1]))
	}
	R0 <- numeric()
	r  <- numeric()
	i0 <- numeric()
	
	for(i in seq_along(unique(ts$mc))){
		print(i)
		tmp   <- estim.R0.impl.one.mc(i, t.max.fit, ts, gi.bck.mean.mc)
		R0[i] <- tmp$R0
		r[i]  <- tmp$r
		i0[i] <- tmp$i0
	}
	
	R0.avg <- mean(R0)
	r.avg  <- mean(r)
	i0.avg <- mean(i0)
	
	# Calculate average exponential growth in incidence:
	ts.avg$ig.avg <- exp(r.avg * ts.avg$time + i0.avg) 
	ts.avg$ig.avg[ts.avg$time>t.max.fit] <- NA
	
	# plot
	if(do.plot){
		g <- ggplot(ts.avg) + geom_line(aes(x=time,y=inc.m)) #+ geom_point(aes(x=time,y=inc.m)) 
		g <- g + scale_y_log10()
		g <- g + geom_line(aes(x=time,y=ig.avg), colour='blue', size=2, alpha=0.5)
		g <- g + geom_ribbon(aes(x=time, ymin=inc.lo, ymax=inc.hi), alpha=0.3)
		g <- g + geom_vline(xintercept = t.max.fit, linetype=2)
		g <- g + ggtitle(paste('Mean incidence time series ; Implied R0 SIR =', round(R0.avg,3))) + ylab('Mean Incidence')
		plot(g)
		g <- ggplot(data = data.frame(R0=R0)) + geom_histogram(aes(x=R0),binwidth = 0.2)
		g <- g + ggtitle(paste('Mean R0 = ',round(R0.avg,3))) + xlab('')
		g <- g + geom_vline(aes(xintercept=R0.avg), size=3, colour='red')
		plot(g)
	}
	return(R0.avg)
}



identify.fizzle <- function(pop.all.mc){
	# Return the MC iterations that were fizzles
	x <- ddply(pop.all.mc,c('mc'),summarize, fz = sum(is_recovered), n=length(id_indiv))
	
	x$r <- x$fz / x$n
	x$fizzle <- FALSE
	x$fizzle[x$r<0.03] <- TRUE
	
	fizz <- x$fizzle
	names(fizz) <- x$mc
	return(fizz)
}


# Calculate the proportion of fizzled simulations
proportion.fizzles <- function(pop) {
	a <- ddply(pop,c('mc'),summarize, 
			   n = length(unique(id_indiv)),
			   m = sum(is_recovered))
	# Final size as a proportion
	a$r <- a$m / a$n
	# Identify fizzle:
	a$is.fizzle <- FALSE
	thresh <- 0.11
	a$is.fizzle[a$r < thresh] <- TRUE
	# Proportion of fizzles:
	return(list(p=sum(a$is.fizzle)/nrow(a),
				n=nrow(a))
	)
}

plot.prop.fizzles <- function(pop){
	p <- proportion.fizzles(pop)
	barplot(p$p,ylim = c(0,1), 
			main=paste('Proportion of fizzles =',
					   round(p$p,3),
					   '(n.MC =',p$n,')'))
	abline(h=1)
}


# Calculate vaccine efficacy on various outcomes
vax.efficacy <- function(pop, outcome){
	
	u   <- subset(pop, is_vaccinated==0)
	v   <- subset(pop, is_vaccinated==1)
	
	# Total number:
	n.u <- nrow(u)
	n.v <- nrow(v)
	
	# Count individual with outcome:
	cnt.u <- sum(u[[outcome]])
	cnt.v <- sum(v[[outcome]])
	
	# proportions:
	p.u <- cnt.u/n.u
	p.v <- cnt.v/n.v
	
	efficacy <- 1 - p.v/p.u
	return(efficacy)
}

plot.vax.efficacy <- function(pop) {
	outcome <- c('is_recovered', 'was_symptomatic', 'was_hosp', 'is_alive')
	eff     <- sapply(outcome, vax.efficacy, pop=pop)	
	df      <- data.frame(outcome, efficacy=eff)
	
	g <- ggplot(df)+geom_bar(aes(x=outcome, y=efficacy), stat='identity')
	g <- g + ggtitle('Vaccine efficacy') #+ coord_cartesian(ylim=c(0,1))
	return(g)
}



# Create a synthetic age distribution for adults.
synthetic_age_adult <- function(age.adult){
	age.thres <- 60
	age.max <- max(age.adult)
	idx <- which(age.adult > age.thres)
	# Decline of older adults:
	rel.prop <- (age.max - age.adult[idx])/(age.max - age.thres)
	
	p.adult <- c( rep(1.0, idx[1]-1) , rel.prop)
	p.adult <- p.adult/sum(p.adult)
	return(p.adult)
}



merge.pop.mc <- function(res.list, n.cpu, 
						 doparallel = FALSE, 
						 select.mc = NA) {
	# Merge all MC iterations into one single data frame for population.
	
	snwrap <- function(i){
		res    <- res.list[[i]]
		pop    <- as.data.frame(res[['world_final']])
		pop$mc <- i
		return(pop)
	}
	t00<- as.numeric(Sys.time())
	seqmc <- seq_along(res.list)
	if(!is.na(select.mc[1])) seqmc <- select.mc
	
	sfInit(parallel = doparallel, cpu = n.cpu) 
	pop.list <- list()
	print('Merging populations...')
	pop.list <- sfSapply(seqmc, snwrap, simplify= FALSE)
	sfStop()
	
	pop.all <- do.call('rbind.data.frame',pop.list)
	t01  <- as.numeric(Sys.time())
	msgt <- round((t01-t00)/60,2)
	print(paste('... populations merged in',msgt,'minutes.'))
	return(pop.all)
}


merge.ts.mc <- function(res.list, n.cpu=2, 
						is.contact=FALSE, is.sp=FALSE,
						doparallel = FALSE){
	# Merge all MC iterations into one single data frame for population.
	tt1 <- as.numeric(Sys.time())
	print('Merging time series...')
	if(is.contact) print('(contacts)')
	if(is.sp) print('(social places)')
	
	snwrap <- function(i){
		print(paste(i,'/',length(res.list)))
		res   <- res.list[[i]]
		if(is.contact)  ts    <- as.data.frame(res[['track_n_contacts']])
		if(!is.contact) ts    <- as.data.frame(res[['time_series']])
		if(is.sp)       ts    <- as.data.frame(res[['time_series_sp']])
		ts$mc <- i
		return(ts)
	}
	### ***WARNING***
	### PARALLEL TURNED OFF BECAUSE DOES NOT WORK (FIX THIS IF TIME...)
	sfInit(parallel = FALSE, cpu = n.cpu)
	ts.list <- list()
	ts.list <- sfSapply(seq_along(res.list), snwrap, simplify= FALSE)
	sfStop()
	ts.all <-  do.call('rbind.data.frame',ts.list)
	tt2 <- as.numeric(Sys.time())
	print(paste('... time series merged in',round((tt2-tt1)/60,1),'minutes'))
	return(ts.all)
}


average.age.contact <- function(res.list){
	
	A.list <- list()
	print('Averaging age-contact matrices...')
	for(l in seq_along(res.list)){
		print(paste(l,'/',length(res.list)))
		x <- res.list[[l]]$wiw_ages
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
		A.list[[l]] <- log(1+A)
	}
	
	# make sure all matrices (one for each MC)
	# are of the same size. If not, crop the largest ones.
	mxsz <- 1E6
	for(l in seq_along(A.list)){
		mxsz <- min(mxsz, nrow(A.list[[l]]))
	}
	for(l in seq_along(A.list)){
		A.list[[l]] <- A.list[[l]][1:mxsz,1:mxsz]
	}
	
	# Mean calculation:
	A.mean <- apply(simplify2array(A.list), 1:2, mean)
	print('... matrices averaged.')
	return(list(A.mean,m.max))
}

### Difference between 2 vectors 
### not necessarily the same size.
diff.flex <- function(x,y){
	nx <- length(x)
	ny <- length(y)
	res <- numeric(max(nx,ny))
	
	if(nx >ny){
		res[1:ny] <- x[1:ny] - y
		res[(ny+1):nx] <- x[(ny+1):nx]
	}
	else if(nx<ny){
		res[1:nx] <- x - y[1:nx]
		res[(nx+1):ny] <- -y[(nx+1):ny]
	}
	else if(nx==ny){
		res = x-y
	}
	return(res)
}

sp.to.string <- function(i) {
	# WARNING: order matters!
	# MUST be same order as enum SPType definition (C++).
	if (i==0) res = 'household';
	if (i==1) res = 'workplace'; 
	if (i==2) res = 'school';
	if (i==3) res = 'other';
	if (i==4) res = 'hospital';
	if (i==5) res = 'pubTransp';
	# // [add here newly defined SPs...]
	# // if(i==6) res = SP_xxx;
	return(res)
}



###   =====   P O P U L A T I O N    P L O T S =====

plot.binomial.regression <- function(dat, xvar, binomial_response, title, split.mc=FALSE) {
	n <- nrow(dat)
	
	g <- ggplot(dat) 
	if(n<1000) g <- g + geom_point(aes_string(x=xvar,y=binomial_response), alpha=0.3) 
	if(!split.mc) g <- g + geom_smooth(aes_string(x=xvar, y=binomial_response), 
									   method = "glm", 
									   method.args = list(family = "binomial"), 
									   colour='red3',size=2,se = F)
	if(split.mc) {
		dat$mc <- as.factor(dat$mc)
		g <- g + geom_smooth(aes_string(x=xvar, y=binomial_response, colour='mc'), 
							 method = "glm", 
							 method.args = list(family = "binomial"), 
							 alpha = 0.6,
							 size=2,se = F)
	}
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

# DELETE, OBSOLETE???:
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


plot.age.contact.matrix.avg <- function(res.list) {
	
	# Desired contact assortativity:
	y <- res.list[[1]]$contactAssort
	D <- matrix(unlist(res.list[[1]]$contactAssort), ncol = length(y))
	
	# Effective contacts from simulation:
	aac <- average.age.contact(res.list)
	A.plot <- aac[[1]]
	m.max <- aac[[2]]
	
	# -- PLOTS --
	
	par(mfrow = c(1,2), cex.lab=1.7,cex.axis=1.7)
	
	na <- ncol(A.plot)
	D.plot <- D[1:(na+1),1:(na+1)]
	
	image(x=0:na, y=0:na, z=D.plot, 
		  col = topo.colors(12), 
		  main='Input contact assortativity',
		  xlab = 'age', ylab='age', las=1)
	abline(a=0,b=1,lty=2); grid()
	
	image(A.plot, x = 1:nrow(A.plot), y=1:nrow(A.plot), 
		  zlim = c(0,max(A.plot)), 
		  ylab = 'infector\'s age',
		  xlab = 'infectee\'s age',
		  las = 1,
		  main = 'Simulated effective contact age matrix\n(averaged across all MC iterations)',
		  col  = topo.colors(12))
	abline(a=0,b=1,lty=2); grid()
}



plot.n.contacts <- function(nc){
	df0 <- data.frame(time=nc$time, uid=nc$uid, n=nc$nContacts, mc=nc$mc)
	df0$timeround <- ceiling(df0$time)
	
	df.ts <- ddply(df0, c('timeround','mc'), summarize, tot.contacts = sum(n))
	
	gts <- ggplot(df.ts, aes(x=timeround, y=tot.contacts, dummy=factor(mc)))
	gts <- gts + geom_step(size=2, alpha=0.5)
	gts <- gts + scale_y_log10() + ggtitle('Total number of contacts\n(no fizzles)')+xlab('time')+ylab('')
	
	# Distribution
	m   <- mean(df$ncontacts)
	ci  <- 0.95
	qlo <- quantile(df$ncontacts, probs=(1-ci)/2)
	qhi <- quantile(df$ncontacts, probs=0.5 + ci/2)
	
	g <- ggplot(df,aes(x=ncontacts)) + geom_histogram(binwidth=1,fill='gold2',colour='gold3') 
	g <- g + ggtitle(paste0('Distribution for the number of contacts (no fizzles)\n (mean = ',round(m,2),' ; ',
							ci*100,'%CI: ',round(qlo,2),' -- ',round(qhi,2),')'))
	g <- g + xlab('Contacts per day, per individual')#+scale_y_log10()
	g <- g + geom_vline(xintercept=m, linetype=1, colour='gold3', size=2)
	
	grid.arrange(gts,g, ncol=2)	
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


plot.age.distrib.mc <- function(pop){
	g <- ggplot(pop)
	g.age <- g + geom_density(aes(x=age, colour=factor(mc)))
	g.age <- g.age + ggtitle('Age distribution')
	return(g.age)
}



plot.symptomatic.fraction <- function(pop){
	pop <- ddply(pop.all.mc, c('mc'),summarize, 
				 nis = sum(was_symptomatic), 
				 ni = sum(is_recovered),
				 n = length(id_indiv))
	
	pop$r <- pop$nis / pop$ni
	pop
	
	p.s <- median(pop$r)
	p.a <- 1-p.s
	
	df <- data.frame(y=c(p.a,p.s), inf_type=c('asymptomatic','symptomatic'))
	
	g <- ggplot(df, aes(x=factor(1),y=y,fill=inf_type))  + geom_bar(width = 1,stat = 'identity')+ coord_polar(theta='y')
	g <- g + xlab('')+ylab('')
	g <- g + ggtitle(paste0('Asymptomatic fraction (median = ',round(p.a,2),')'))
	g <- g + scale_fill_brewer(palette = 'Set2')
	return(g)
}

plot.proportion <- function(pop) {
	
	pop$dead <- 1 - pop$is_alive
	df <- ddply(pop, c('mc'), summarize, 
				death     = sum(dead),
				finalsize = sum(is_recovered),
				sympt     = sum(was_symptomatic),
				hosp      = sum(was_hosp),
				popsize   = length(id_indiv))
	
	df$finalsize.prop <- df$finalsize / df$popsize
	df$sympt.prop     <- df$sympt / df$popsize
	df$hosp.prop      <- df$hosp / df$popsize
	df$death.prop     <- df$death / df$popsize
	df$death.hosp.prop<- df$death / df$death
	
	# get rid of fizzles:
	df <- subset(df, finalsize.prop > 0.01)
	df <- df [,grepl('.prop', names(df))]
	
	z <- gather(df, type, 'p',1:5)
	g <- ggplot(z) + geom_boxplot(aes(x=type, y=p, fill=type)) 
	g <- g + facet_wrap(~type, scales = 'free')
	g <- g + theme(axis.title.x=element_blank(),
				   axis.text.x=element_blank(),
				   axis.ticks.x=element_blank())
	return(g)
}




plot.population <- function(pop, split.mc = TRUE) {
	
	pop.hosp <- subset(pop, was_hosp>0)
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
												title='Survival and Frailty',
												split.mc=split.mc)
	
	g.death.imm.hum <- plot.binomial.regression(dat = pop, 
												xvar='immunity_hum',
												binomial_response = 'is_alive',
												title='Survival and Humoral Immunity',
												split.mc=split.mc)
	
	g.death.imm.cell <- plot.binomial.regression(dat = pop, 
												 xvar='immunity_cell',
												 binomial_response = 'is_alive',
												 title='Survival and Cellular Immunity',
												 split.mc=split.mc)
	
	
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
	
	g.dol.drawn <- g + geom_histogram(aes(x=dol_drawn), fill='orange',colour='orange', size=linew, alpha=alpha, bins=30)
	g.dol.drawn <- g.dol.drawn + geom_vline(xintercept = m_dol_drawn, colour='orange',linetype = 2)
	g.dol.drawn <- g.dol.drawn + ggtitle(paste0("DOL drawn distribution, among infected (mean = ", round(m_dol_drawn,2),")"))
	
	g.doi.drawn <- g + geom_histogram(aes(x=doi_drawn), fill='red',colour='red', size=linew, alpha=alpha, bins=30)
	g.doi.drawn <- g.doi.drawn + geom_vline(xintercept = m_doi_drawn, colour='red',linetype = 2)
	g.doi.drawn <- g.doi.drawn + ggtitle(paste0("DOI drawn distribution, among infected (mean = ", round(m_doi_drawn,2),")"))
	
	
	m_dobh_drawn = mean(pop.hosp$dobh_drawn)
	m_doh_drawn  = mean(pop.hosp$doh_drawn)
	
	g.hosp <- ggplot(pop.hosp)
	g.dobh.drawn <- g.hosp + geom_histogram(aes(x=dobh_drawn), fill='purple',colour='purple', size=linew, alpha=alpha,bins=30)
	g.dobh.drawn <- g.dobh.drawn + geom_vline(xintercept = m_dobh_drawn, colour='purple',linetype = 2)
	g.dobh.drawn <- g.dobh.drawn + ggtitle(paste0("DOBH drawn distribution, among hospitalized (mean = ", round(m_dobh_drawn,3),")"))
	
	g.doh.drawn <- g.hosp + geom_histogram(aes(x=doh_drawn), fill='green4',colour='green4', size=linew, alpha=alpha, bins=30)
	g.doh.drawn <- g.doh.drawn + geom_vline(xintercept = m_doh_drawn, colour='green4',linetype = 2)
	g.doh.drawn <- g.doh.drawn + ggtitle(paste0("DOH drawn distribution, among hospitalized (mean = ", round(m_doh_drawn,3),")"))
	
	
	### ==== Hospitalizations ==== 
	
	g.frail.hosp <- plot.binomial.regression(dat = pop, xvar = 'frailty',
											 binomial_response = 'hosp',
											 title = 'Hospitalization and frailty',
											 split.mc=split.mc)
	
	g.imm.hum.hosp <- plot.binomial.regression(dat = pop, xvar = 'immunity_hum',
											   binomial_response = 'hosp',
											   title = 'Hospitalization and humoral immunity',
											   split.mc=split.mc)
	
	g.imm.cell.hosp <- plot.binomial.regression(dat = pop, xvar = 'immunity_cell',
												binomial_response = 'hosp',
												title = 'Hospitalization and cellular immunity',
												split.mc=split.mc)
	
	g.age.hosp <- plot.binomial.regression(dat = pop, xvar = 'age',
										   binomial_response = 'hosp',
										   title = 'Hospitalization and age',
										   split.mc=split.mc)
	
	
	
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
	
	
	# Summarize data
	pop$is_dead <- 1-pop$is_alive
	tmp <- ddply(pop,c('age'),summarize,
				 nh = sum(was_hosp),
				 nd = sum(is_dead),
				 n = length(id_indiv))
	tmp$prop.hosp      <- tmp$nh/tmp$n
	tmp$prop.dead      <- tmp$nd/tmp$n
	tmp$prop.dead.hosp <- tmp$nd/tmp$nh
	
	g.age.hosp.raw <- ggplot(tmp)+geom_bar(aes(x=age, y=prop.hosp),stat = 'identity')
	g.age.hosp.raw <- g.age.hosp.raw + ggtitle('Hospitalizations by age')
	g.age.hosp.raw <- g.age.hosp.raw + xlab('age')+ylab('proportion')
	
	g.age.dead.raw <- ggplot(tmp)+geom_bar(aes(x=age, y=prop.dead),stat = 'identity')
	g.age.dead.raw <- g.age.dead.raw + ggtitle('Whole Population Death ratio by age')
	g.age.dead.raw <- g.age.dead.raw + xlab('age')+ylab('proportion')
	
	g.age.dead.hosp.raw <- ggplot(tmp)+geom_bar(aes(x=age, y=prop.dead),stat = 'identity')
	g.age.dead.hosp.raw <- g.age.dead.hosp.raw + ggtitle('Hospitalized Deaths ratio by age')
	g.age.dead.hosp.raw <- g.age.dead.hosp.raw + xlab('age')+ylab('proportion')
	
	### ==== Secondary cases distribution ==== 
	R0    <- mean(pop$n_secondary_cases)
	nsmax <- max(pop$n_secondary_cases)
	g.R <- ggplot(pop) + geom_histogram(aes(n_secondary_cases), 
										fill = 'red3',
										colour = 'red4',
										breaks = seq(0,nsmax+1,by=1))
	g.R <- g.R + scale_x_continuous(breaks =  seq(0,nsmax+1,by=1))
	g.R <- g.R + geom_vline(aes(xintercept=R0), linetype = 2, size = 2)
	g.R <- g.R + ggtitle(paste0("Secondary cases distribution (R0=",
								round(R0,2),")"))
	
	# R0 symptomatic or not
	pop.transm.sum <- ddply(pop,"was_symptomatic",
							summarize,
							m = mean(n_secondary_cases),
							q.lo = quantile(n_secondary_cases,probs = 0.5-0.80/2),
							q.hi = quantile(n_secondary_cases,probs = 0.5+0.80/2))
	
	pop.transm.sum
	
	g.R.symptom <- ggplot(pop) +geom_histogram(aes(x=n_secondary_cases, ..density..,
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
	
	pop.transm.sum2 <- ddply(pop,"is_treated",
							 summarize,
							 m = mean(n_secondary_cases),
							 q.lo = quantile(n_secondary_cases,probs = 0.5-0.80/2),
							 q.hi = quantile(n_secondary_cases,probs = 0.5+0.80/2))
	
	pop.transm.sum2
	
	g.R.treat <- ggplot(pop) +geom_histogram(aes(x=n_secondary_cases, ..density..,
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
											 bins=30)
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
	
	g.sympt.frac <- plot.symptomatic.fraction(pop)
	
	g.sympt.imm.hum.dist <- plot.density.categ(dat = pop, 
											   xvar ='immunity_hum',
											   categ = 'was_symptomatic',
											   title ='Humoral Immunity distribution\n(by symptomatic status)')
	
	g.sympt.imm.cell.dist <- plot.density.categ(dat = pop, 
												xvar ='immunity_cell',
												categ = 'was_symptomatic',
												title ='Cellular Immunity distribution\n(by symptomatic status)')
	
	
	sympt.vax <- ddply(pop.inf,c("is_vaccinated","was_symptomatic"), summarize, n=length(id_indiv))
	
	# Vaccine efficacy:
	# u   <- subset(sympt.vax, is_vaccinated==0)
	# u.s <- subset(sympt.vax, is_vaccinated==0 & was_symptomatic==1)
	# prop.sympt.no.vax <- u.s$n / sum(u$n)
	# 
	# v   <- subset(sympt.vax, is_vaccinated==1)
	# v.s <- subset(sympt.vax, is_vaccinated==1 & was_symptomatic==1)
	# prop.sympt.vax <- v.s$n / sum(v$n)
	# 
	# vax.eff <- 1 - prop.sympt.vax/prop.sympt.no.vax
	
	g.sympt.vax <- ggplot(sympt.vax)+geom_bar(aes(x=factor(is_vaccinated), y=n, fill=factor(was_symptomatic)),
											  stat='identity', position='dodge')
	g.sympt.vax <- g.sympt.vax + ggtitle('Symptomatic infection and vaccination')
	
	g.vax.eff <- plot.vax.efficacy(pop)
	
	# === Schedules ===
	
	g.sched.R <- plot.sched(pop)
	
	### ==== Final ====
	
	grid.arrange(g.age, 
				 g.age.imm.hum,
				 g.age.imm.cell,
				 g.age.fra,
				 g.age.death.dist)
	
	grid.arrange(g.death.frailty,
				 g.death.imm.hum,
				 g.death.imm.cell,
				 g.death.frailty.dist,
				 g.death.imm.hum.dist,
				 g.death.imm.cell.dist)
	
	grid.arrange(g.age.hosp.raw,
				 g.age.dead.raw,
				 g.age.dead.hosp.raw,
				 g.imm.hum.hosp,
				 g.imm.cell.hosp,
				 g.frail.hosp,
				 g.age.hosp, 
				 g.treat.hosp,
				 g.vax.hosp)
	
	grid.arrange(g.sympt.frac,
				 g.sympt.vax,
				 g.vax.eff,
				 g.sympt.imm.hum.dist,
				 g.sympt.imm.cell.dist,
				 g.sched.R)
	
	try(grid.arrange( g.dol.drawn, 
					  g.doi.drawn,
					  g.dobh.drawn,
					  g.doh.drawn), silent = TRUE)
	
	grid.arrange( g.R,
				  g.R.symptom,
				  g.R.treat,
				  g.gibck,
				  g.gibck.sympt)
	
}



###   =====   S O C I A L   P L A C E S   P L O T S =====

plot.world <- function(x) {
	stopifnot(length(x)>=1)
	sp      <- as.data.frame(x[[1]]$census_sp)
	stopifnot(nrow(sp)>1)
	sp$type <- sapply(sp$sp_type,FUN = sp.to.string) 
	sp.cnt  <- ddply(sp,c('type'),summarize, cnt=length(sp_id))
	
	sp.cnt$target <- NA
	sp.cnt$target[sp.cnt$type=='household'] <- world.prm[['n_hh']]
	sp.cnt$target[sp.cnt$type=='workplace'] <- world.prm[['n_wrk']]
	sp.cnt$target[sp.cnt$type=='other']     <- world.prm[['n_other']]
	sp.cnt$target[sp.cnt$type=='school']    <- world.prm[['n_school']]
	sp.cnt$target[sp.cnt$type=='pubTransp'] <- world.prm[['n_pubt']]
	sp.cnt$target[sp.cnt$type=='hospital']  <- world.prm[['n_hosp']]
	sp.cnt$prop.filled <- sp.cnt$cnt/sp.cnt$target
	sp.cnt
	
	g <- ggplot(sp.cnt) + geom_bar(aes(x=type, y=prop.filled), stat='identity')
	g <- g + ggtitle('Proportion filled w.r.t. target sizes')
	
	sp.dsz  <- ddply(sp,c('type','nlinked'), summarize, n=length(nlinked))
	
	q <- join(sp.dsz, sp.cnt,by='type')
	q$prop <- q$n/q$cnt
	q$source <- 'simulated'
	
	wp <- data.frame(type = 'household', nlinked=world.prm[['hh_size']], prop=world.prm[['hh_size_proba']])
	wp <- rbind(wp, 
				data.frame(type = 'workplace', nlinked=world.prm[['wrk_size']], prop=world.prm[['wrk_size_proba']]))
	wp <- rbind(wp, 
				data.frame(type = 'school', nlinked=world.prm[['school_size']], prop=world.prm[['school_size_proba']]))
	wp <- rbind(wp, 
				data.frame(type = 'pubTransp', nlinked=world.prm[['pubt_size']], prop=world.prm[['pubt_size_proba']]))
	wp <- rbind(wp, 
				data.frame(type = 'other', nlinked=world.prm[['other_size']], prop=world.prm[['other_size_proba']]))
	
	wp$type <- as.character(wp$type)
	wp$source <- 'target'
	
	df <- rbind(q[,names(wp)],wp)
	df$lt <- 1
	df$lt[df$source=='target'] <- 2
	g2 <- ggplot(df,aes(x=nlinked, y=prop, colour=source, shape=source)) + geom_point() + geom_line(aes(linetype=factor(lt))) + facet_wrap(~type, scale='free')
	g2 <- g2 + ggtitle('Size Distribution: simulated v.s. target')
	
	plot(g)
	plot(g2)
}


###   =====   T I M E  S E R I E S     P L O T S =====

plot.epi.timeseries <- function(ts){
	ts$death_incidence <- c(ts$nD[1],diff(ts$nD))
	
	idx <- !grepl('time',names(ts))
	idx <- idx * !grepl('mc',names(ts))
	idx <- idx * c(1:length(idx))
	idx <- idx[idx>0]
	tsl <- gather(data = ts, key = "type",value = "n",idx)
	tsl$day <- round(tsl$time)
	
	D <- ddply(tsl,c('day','type'),summarise, 
			   m = median(n),
			   q.lo = quantile(n,probs = 0.025),
			   q.hi = quantile(n,probs = 0.975))
	
	plot.ts <- function(D, do.log, title='') {
		g <- ggplot(D) + geom_line(aes(x=day, y=m+0.1, colour=type),size=1)
		g <- g + geom_ribbon(aes(x=day,ymin=q.lo,ymax=q.hi,fill=type),alpha=0.2)
		if(do.log) g <- g + scale_y_log10()
		g <- g + ylab('') + ggtitle(title)
		return(g)
	}
	
	D.SR <- subset(D, type %in% c('nS','nR'))
	D.E <- subset(D, type %in% c('nE'))
	D.I <- subset(D, type %in% c('nIa','nIs','prevalence'))
	
	D.i.vt <- subset(D, type %in% c('incidence','n_treated','n_vaccinated'))
	D.d.vt <- subset(D, type %in% c('death_incidence','n_treated','n_vaccinated'))
	
	grid.arrange(
		plot.ts(D.SR, do.log=F, title='Susceptible & Recovered'),
		plot.ts(D.E, do.log=T, title='Exposed'),
		plot.ts(D.I, do.log=F,title = 'Prevalence'),
		plot.ts(D.I, do.log=T,title = 'Prevalence (log-scale)'),
		plot.ts(D.i.vt, do.log=T, title = 'Incidence, vaccination & treatment'),
		plot.ts(D.d.vt, do.log=T, title = 'Deaths, vaccination & treatment')
	)
}


plot.epi.timeseries.comp <- function(ts){
	ts$death_incidence <- c(ts$nD[1],diff(ts$nD))
	
	idx <- !grepl('time',names(ts))
	idx <- idx * !grepl('mc',names(ts))
	idx <- idx * !grepl('scen',names(ts))
	idx <- idx * c(1:length(idx))
	idx <- idx[idx>0]
	tsl <- gather(data = ts, key = "type",value = "n",idx)
	tsl$day <- round(tsl$time)
	
	D <- ddply(tsl,c('day','type','scen'),summarise, 
			   m = median(n),
			   q.lo = quantile(n,probs = 0.025),
			   q.hi = quantile(n,probs = 0.975))
	
	plot.ts <- function(D, do.log, title='') {
		g <- ggplot(D) + geom_line(aes(x=day, y=m+0.1, colour=scen),
								   size=1)
		g <- g + geom_ribbon(aes(x=day,ymin=q.lo,ymax=q.hi, fill=scen),alpha=0.2)
		if(do.log) g <- g + scale_y_log10()
		g <- g + ylab('') + ggtitle(title)
		g <- g + facet_wrap(~type,scales = 'free')
		return(g)
	}
	
	D.SR <- subset(D, type %in% c('nS','nR'))
	D.E <- subset(D, type %in% c('nE'))
	D.I <- subset(D, type %in% c('nIa','nIs','prevalence'))
	D.cuminc <- subset(D, type %in% c('cuminc'))
	
	D.i.vt <- subset(D, type %in% c('incidence','nD'))
	D.d.vt <- subset(D, type %in% c('n_treated','n_vaccinated'))
	
	grid.arrange(
		plot.ts(D.SR, do.log=F, title='Susceptible & Recovered'),
		plot.ts(D.E, do.log=T, title='Exposed'),
		plot.ts(D.I, do.log=F,title = 'Prevalence'),
		plot.ts(D.I, do.log=T,title = 'Prevalence (log-scale)'),
		plot.ts(D.i.vt, do.log=F, title = 'Incidence & Cumul Deaths'),
		plot.ts(D.d.vt, do.log=F, title = 'Vaccination & treatment')
	)
}


plot.ts.sp <- function(tssp){
	
	tssp$timeround <- floor(tssp$time)
	df0 <- ddply(tssp,c('timeround','type','mc'),summarize, n=sum(nE))
	
	df <- ddply(df0,c('timeround','type'),summarize, 
				m=median(n), 
				q.lo=quantile(n,probs=0.025),
				q.hi=quantile(n,probs=0.975))
	
	g <- ggplot(df,aes(x=timeround))
	g <- g +geom_line(aes(colour=type, y=m),size=1) #+ geom_point(aes(colour=type, y=m))
	g <- g + geom_ribbon(aes(ymin=q.lo, ymax=q.hi, fill=type),alpha=0.2)
	g <- g + ggtitle("Exposed individuals, by SP types") + xlab('day')
	g <- g + scale_y_log10()
	
	gf <- g +facet_wrap(~type)
	grid.arrange(g,gf)
}


sp_type_string <- function(sp_type_int) {
	if(sp_type_int==0) return('Household')
	if(sp_type_int==1) return('Workplace')
	if(sp_type_int==2) return('School')
	if(sp_type_int==3) return('Other')
	if(sp_type_int==4) return('Hospital')
	if(sp_type_int==5) return('PubTransp')
	return(NA)
}


read_input_sp_size_distribution <- function(world.prm) {
	# Read target distributions and reformat in data frame:
	hh.x     <- world.prm$hh_size
	hh.y     <- world.prm$hh_size_proba
	wrk.x    <- world.prm$wrk_size
	wrk.y    <- world.prm$wrk_size_proba
	school.x <- world.prm$school_size
	school.y <- world.prm$school_size_proba
	pubt.x   <- world.prm$pubt_size
	pubt.y   <- world.prm$pubt_size_proba
	other.x  <- world.prm$other_size
	other.y  <- world.prm$other_size_proba
	
	df <- data.frame(size=NULL, freq=NULL, sp_type_string=NULL)	
	
	for(i in 1:length(hh.x)){
		df.tmp <- data.frame(size=hh.x[[i]], freq=hh.y[[i]], sp_type_string='Household')	
		df <- rbind(df,df.tmp)
	}
	for(i in 1:length(wrk.x)){
		df.tmp <- data.frame(size=wrk.x[[i]], freq=wrk.y[[i]], sp_type_string='Workplace')	
		df <- rbind(df,df.tmp)
	}
	for(i in 1:length(school.x)){	
		df.tmp <- data.frame(size=school.x[[i]], freq=school.y[[i]], sp_type_string='School')	
		df <- rbind(df,df.tmp)
	}
	for(i in 1:length(pubt.x)){
		df.tmp <- data.frame(size=pubt.x[[i]], freq=pubt.y[[i]], sp_type_string='PubTransp')	
		df <- rbind(df,df.tmp)
	}
	for(i in 1:length(other.x)){
		df.tmp <- data.frame(size=other.x[[i]], freq=other.y[[i]], sp_type_string='Other')	
		df <- rbind(df,df.tmp)
	}
	df$sp_type_string <- as.character(df$sp_type_string)
	return(df)
}


plot.sp.sz.distrib <- function(res.list.0, world.prm){
	
	par(mfrow=c(2,2))
	
	# Read target data:
	df <- read_input_sp_size_distribution(world.prm)
	# Retrieve all SocialPlaces types:
	sptype.vec <- seq_along(res.list.0[[1]][['sp_size_distrib']])
	
	
	for(k in seq_along(sptype.vec)) {
		
		sptype <- sp_type_string(k-1)
		data.df <- subset(df, sp_type_string==sptype)
		
		if(sptype != 'Hospital' & sptype != 'Other'){
			
			# Retrieve social place sizes of the kth type, 
			# for the mth MC iteration:
			p <- list()
			cnt <- 1
			for( m in idx.mc.no.fizz){
				sp.sz.d <- res.list.0[[m]][['sp_size_distrib']]
				h <- hist(x = sp.sz.d[[k]], plot = F, breaks = c(0,data.df$size))
				p[[cnt]] <- h$counts / sum(h$counts)
				cnt <- cnt + 1
			}
			# Summary stats across all MC iterations:
			pm <- matrix(unlist(p), ncol=cnt-1)
			p.avg <- apply(X = pm,MARGIN = 1,FUN = mean)
			p.max <- apply(X = pm,MARGIN = 1,FUN = max)
			p.min <- apply(X = pm,MARGIN = 1,FUN = min)
			
			# Plots:
			plot(x = data.df$size, y = data.df$freq, 
				 ylim = range(data.df$freq, p.avg),
				 col='red',lwd=1, cex=2, typ='o', lty=2, las=1, 
				 main = sptype, xlab = 'size', ylab='proportion')
			grid()
			lines(x = data.df$size, y = p.avg, typ='o', pch=16, cex=1.7)
			lines(x = data.df$size, y = p.max, typ='l', lty=3)
			lines(x = data.df$size, y = p.min, typ='l', lty=3)
			legend(x='topright', legend = c('target','simulation'), col=c('red','black'), pch=c(1,16), lty=c(2,1))
		}
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


plot.ts.comp <- function(df,var.name) {
	df$day <- round(df$time,0)
	df$var <- df[,var.name]
	df2 <- ddply(df,c('day'),summarize, 
				 m = mean(var),
				 md = median(var),
				 q.lo = quantile(var,probs = 0.025),
				 q.hi = quantile(var,probs = 0.975))
	
	g <- ggplot(df2) + geom_line(aes(x=day,y=m),size=1) + geom_line(aes(x=day,y=md),linetype=2)  
	g <- g + geom_ribbon(aes(x=day,ymin=q.lo,ymax=q.hi),alpha=0.2)
	g <- g + ggtitle(var.name)
	g
}

plot.ts.comp.all <- function(df){
	grid.arrange(
		plot.ts.comp(df,'dinc'),
		plot.ts.comp(df,'dIs'),
		plot.ts.comp(df,'dIa')
	)
}





