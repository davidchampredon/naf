### 
###   COMPARE SIMULATIONS (BASELINE vs INTERVENTION) ALREADY RUN
###

ct0 <- as.numeric(Sys.time())
library(tidyr)

### ==== Load simulation results ====

print('Loading simulation results ...')
load('mc-simul.RData')
t1 <- as.numeric(Sys.time())
print(paste('... simulation results loaded in',round((t1-ct0)/60,1),'minutes.'))

if(is.na(res.list)){
	print('NOTHING TO COMPARE... EXITING.')
	quit()
}

n.mc  <- length(res.list)
print(paste('Number of MC iterations:',n.mc))

source('analysis_tools.R')  
save.plot.to.file <- TRUE
max.cpu <- 2
dir.results <- '../results/'

### ==== Merge all MC iterations ====

n.cpu <- min(parallel::detectCores(), max.cpu)
ts    <- merge.ts.mc(res.list,   n.cpu = n.cpu)
ts0   <- merge.ts.mc(res.list.0, n.cpu = n.cpu)

ts0$scen <- 'baseline'
ts$scen  <- 'interv'
u        <- rbind.data.frame(ts0,ts)


tmp     <- list()
dprev   <- numeric(n.mc)
dtreat  <- numeric(n.mc)
dD      <- numeric(n.mc)
dcuminc <- numeric(n.mc)

# Loop to calculate differences
# b/w baseline and intervention:
for(i in seq_along(res.list.0)){
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

dfall <- do.call('rbind.data.frame',tmp)


#### ====== Plots =====

plot.hist <- function(x) {
  hist(x, xlim=range(0,x),
       col='lightgrey',
       #yaxt='n',
       ylab='',xlab='', las=1,
       main=deparse(substitute(x)))
  abline(v=mean(x),lty=1, lwd=4,col='red')
  abline(v=0,lty=2, lwd=2,col='black')
}

pdf(paste0(dir.results,'plot-compare.pdf'), width = 15, height = 10)

plot.epi.timeseries.comp(u)
plot.ts.comp.all(dfall)

par(mfrow=c(1,3))
plot.hist(dtreat)
plot.hist(dD)
plot.hist(dcuminc)

par(mfrow=c(1,1))
CI <- 0.95
df <- data.frame(dtreat,dD,dcuminc)
df0 <- gather(df)
df <- ddply(df0,'key',summarize, 
            m = mean(value), 
            md = median(value),
            q.lo = quantile(value, probs = 0.5-CI/2),
            q.hi = quantile(value, probs = 0.5+CI/2))

g <- ggplot(df) + geom_pointrange(aes(x=key, y=md, ymin=q.lo, ymax=q.hi))
g <- g + geom_point(aes(x=key, y=m),shape=2)
plot(g)
h <- ggplot(df0)+geom_density(aes(x=value),fill='blue')+facet_wrap(~key, scales='free')
plot(h)

dev.off()
t2 <- as.numeric(Sys.time())
print(paste('Full comparison completed in',round((t2-ct0)/60,2),'minutes.'))

