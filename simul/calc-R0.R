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

# Merge all simulations:

n.cpu <- parallel::detectCores() - 1
if(simul.prm[['baseline_only']])  res.select <- res.list.0
if(!simul.prm[['baseline_only']]) res.select <- res.list

pop.all.mc   <- merge.pop.mc(res.select,
							 n.cpu = n.cpu, 
							 doparallel = TRUE)

# METHOD 1: Calculate R0 from realized number of secondary infections:
R0 <- calc.R0(pop.all.mc, time.init = 2)

# METHOD 2: Calculate R0 implied from SIR:
R0.SIR <- calc.R0.SIR(pop.all.mc = pop.all.mc, t.max.fit = 10, do.plot = T)

yy <- numeric()
yy.sir <- numeric()
xx <- c(seq(1.5,3,by=0.25), 4:30)

for(i in seq_along(xx)) {
	yy[i]     <- calc.R0(pop.all.mc, time.init = xx[i])
	yy.sir[i] <- calc.R0.SIR(pop.all.mc = pop.all.mc, t.max.fit = xx[i])
}
plot(xx,yy, typ='o', ylim=c(0,max(yy,na.rm = T)))
lines(xx,yy.sir, typ='o', ylim=c(0,max(yy.sir,na.rm = T)), col='red')
abline(h=1)