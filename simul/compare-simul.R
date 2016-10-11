### 
###   COMPARE SIMULATIONS (BASELINE vs INTERVENTION) ALREADY RUN
###

t0 <- as.numeric(Sys.time())

### ==== Load simulation results ====

print('Loading simulation results ...')
load('mc-simul.RData')
t1 <- as.numeric(Sys.time())
print(paste('... simulation results loaded in',round((t1-t0)/60,1),'minutes.'))

source('analysis_tools.R')  
save.plot.to.file <- TRUE


### ==== Merge all MC iterations ====

n.cpu <- min(parallel::detectCores(), 16)
ts    <- merge.ts.mc(res.list,   n.cpu = n.cpu)
ts0   <- merge.ts.mc(res.list.0, n.cpu = n.cpu)

ts0$scen <- 'baseline'
ts$scen  <- 'interv'
u <- rbind.data.frame(ts0,ts)
plot.epi.timeseries.comp(u)


n.mc    <- length(res.list)
print(paste('Number of MC iterations:',n.mc))
tmp     <- list()
dprev   <- numeric(n.mc)
dtreat  <- numeric(n.mc)
dD      <- numeric(n.mc)
dcuminc <- numeric(n.mc)

for(i in seq_along(res.list.0)){
	print(i)
	
	z0 <- subset(ts0,mc==i)
	z <- subset(ts,mc==i)
	
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
	
	# Cumulative quantities (only last difference relevant):
	dprev[i]   <- (z$prevalence[mn] - z0$prevalence[mn])/pop.size
	dtreat[i]  <- (z$n_treated[mn] - z0$n_treated[mn])/pop.size
	dD[i]      <- (z$nD[mn] - z0$nD[mn])/pop.size
	dcuminc[i] <- (cumsum(z$incidence)[mn] - cumsum(z0$incidence)[mn])/pop.size
	
	tmp[[i]] <- data.frame(time = tmx, 
						   mc   = factor(i),
						   dinc = diff.flex(z$incidence, z0$incidence)/pop.size ,
						   dIa  = diff.flex(z$nIa, z0$nIa)/pop.size,
						   dIs  = diff.flex(z$nIs, z0$nIs)/pop.size
						   )
}

df <- do.call('rbind.data.frame',tmp)


plot.ts.comp.all()


plot.hist <- function(x) {
	hist(x, xlim=range(0,x),
		 col='lightgrey',
		 yaxt='n',ylab='',xlab='',
		 main=deparse(substitute(x)))
	abline(v=mean(x),lty=1, lwd=4,col='red')
	abline(v=0,lty=2, lwd=2,col='black')
}

par(mfrow=c(2,2))
plot.hist(dprev)
plot.hist(dtreat)
plot.hist(dD)
plot.hist(dcuminc)

t2 <- as.numeric(Sys.time())
print(paste('Full comparison completed in',round((t2-t0)/60,1),'minutes.'))

