### 
###   CALCULATE R0 OF SIMULATIONS ALREADY RUN
###

t00 <- as.numeric(Sys.time())

### ==== Load simulation results ====

print('Loading simulation results ...')
load('mc-simul.RData')
t01 <- as.numeric(Sys.time()) ; msgt <- round((t01-t00)/60,2)
print(paste('... simulation results loaded in',msgt,'minutes.'))

source('analysis_tools.R')  

### ==== Merge all MC iterations ====

n.cpu <- parallel::detectCores() - 1
if(simul.prm[['baseline_only']])  res.select <- res.list.0
if(!simul.prm[['baseline_only']]) res.select <- res.list

pop.all.mc   <- merge.pop.mc(res.select,
							 n.cpu = n.cpu, 
							 doparallel = TRUE)

R0 <- calc.R0(pop.all.mc, time.init = 22)

yy <- numeric()
xx <- c(seq(1.5,3,by=0.25), 4:30)
for(i in seq_along(xx)) yy[i] <- calc.R0(pop.all.mc, time.init = xx[i])
plot(xx,yy, typ='o', ylim=c(0,max(yy)))
abline(h=1)


# Merging populations of all MC iterations
# can be very slow, and they very similar
# across MC, so just display the first non-fizzled('select.mc=idx.mc.no.fizz[1]')	
idx.mc         <- unique(pop.all.mc$mc)
fizz.mc        <- identify.fizzle(pop.all.mc)
fizz.mc.idx    <- which(fizz.mc==TRUE)
if(length(fizz.mc.idx)==0) 
	idx.mc.no.fizz <- idx.mc
if(length(fizz.mc.idx)>0) 
	idx.mc.no.fizz <- idx.mc[-fizz.mc.idx]

pop.nofizz <- subset(pop.all.mc, mc %in% idx.mc.no.fizz)

tmp <- subset(pop.nofizz, gi_bck>0)
gi_bck.mean <- mean(tmp$gi_bck)


# Merging time series if faster and more informative,
# so it is done across all MC iterations:
res.no.fizz <- list()
for(i in seq_along(idx.mc.no.fizz)) res.no.fizz[[i]] <- res.list.0[[i]]

ts <- merge.ts.mc(res.no.fizz, n.cpu = n.cpu)
g <- ggplot(ts, aes(x=time,y=incidence+1,colour=factor(mc))) + geom_line()+geom_point()+scale_y_log10()
plot(g)


ts.init <- subset(ts, incidence<30 & nR < 200)


g <- ggplot(ts.init, aes(x=time,y=incidence+1,colour=factor(mc))) + geom_line()+geom_point()+scale_y_log10()
plot(g)


mylm <- function(xx,yy) {
	z <- lm(formula = yy ~ xx)
	return(z$coefficients[2])
}

y <- ddply(ts.init, c('mc'), summarise, r=mylm(time,log(incidence+1)))
y
y$R0 <- 1+y$r*2
y



msg <- paste('R0 = ',round(R0,3))
print(msg); message(msg)
