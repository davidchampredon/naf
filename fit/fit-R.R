###
###   FIT CONTACT RATE TO R0
###
###   Very basic grid search. 
###   Return contact rate values closest to targeted R.

t0 <- as.numeric(Sys.time())

library(snowfall)
library(plyr)

# Load all default parameters:
source('fit_fcts.R')
source('../simul/utils-run.R')
source('../simul/utils-misc.R')
R.library.dir   <- '../Rlibrary/lib'
param.model.dir <- '../param-model/'
data.dir        <- '../data/'
simul.dir       <- '../simul/'

# MC iterations and CPUs sed for this fit:
n.MC  <- 3
n.cpu <- 3	

# Unpack parameters:
PRM <- load.all.parameters(R.library.dir ,
						   param.model.dir,
						   data.dir,
						   simul.dir) 
for(i in seq_along(PRM)){
	assign(x = names(PRM)[i], value = PRM[[i]])
}
stoch_build_world <- simul.prm[['build_world_stoch']]

# Target mean R:
R.target <- 1.8

# Range explored:
cr.mean <- c(1,2,5) #seq(0.5, 6.0, by = 0.25)
cr.cv   <- c(0.5, 0.75)  #c(0.5, 0.75, 1)

n.cpu.cr.mean <- 1

# ==== Simulation Loop ====
# Run the simulations on all scenarios:

wrap_sim <- function (i, cr.mean, cr.cv, 
					  simul.prm, prm, n.MC, n.cpu)
{
	x <- numeric(length = length(cr.cv))
	
	for(j in seq_along(cr.cv)){
		# Overwrite parameter values
		# associated with current scenario:
		prm[['contact_rate_mean']] <- cr.mean[i]
		prm[['contact_rate_CV']]   <- cr.cv[j]

		# Run the simulation for that scenario:
		stoch_build_world <- simul.prm[['build_world_stoch']]
		res.list.0 <- run.simul.fit.R(n.MC, n.cpu, prm, simul.prm, 
									  interv.prm.0, world.prm, sched.prm, 
									  stoch_build_world)	
		
		# merge all MC iterations:
		pop   <- merge.pop.mc(res.list   = res.list.0,
							  n.cpu      = n.cpu, 
							  doparallel = TRUE)
		
		x[j] <- calc.R0.SIR(pop.all.mc = pop, 
							res.list.0 = res.list.0,
							t.max.fit  = 10)  
	}
	return(x)
}

# Parallel computation:
sfInit(parallel = (n.cpu.cr.mean>1) , cpu = n.cpu.cr.mean)
sfLibrary(naf, lib.loc = R.library.dir)
sfLibrary(snowfall)
sfLibrary(plyr)
sfExportAll()
res.tmp <- sfSapply(x = seq_along(cr.mean), fun = wrap_sim, 
					cr.mean=cr.mean, cr.cv=cr.cv,
					simul.prm=simul.prm,
					prm=prm, n.MC=n.MC, 
					n.cpu=n.cpu,
					simplify = FALSE)
sfStop()

res <- matrix(unlist(res.tmp),nrow = length(cr.mean))

save.image(file = 'fit-R.RData')


# ==== Best Fit ====

dist <- abs(res - R.target)
idx <- which(dist == min(dist,na.rm = TRUE), arr.ind = TRUE)
cr.mean.best <- cr.mean[idx[1]]
cr.cv.best   <- cr.cv[idx[2]]

write.csv(x = data.frame(cr.mean.best=cr.mean.best, cr.cv.best=cr.cv.best),
		  file = 'fitted_cr.csv', quote = F, row.names = F)


# ==== PLOTS ====

pdf('fitted_cr.pdf')

plot(1,xlim=range(cr.cv), ylim=range(res,na.rm = TRUE),
	 xlab = 'cv', ylab='R')
for(i in seq_along(cr.mean)){
	lines(x=cr.cv, y=res[i,], typ='o', col=i)
	text(x=cr.cv[1], y=res[i,1],labels = cr.mean[i],pos=4,col=i)
}

plot(1,xlim=range(cr.mean), ylim=range(res,na.rm = TRUE),
	 xlab = 'mean', ylab='R')
for(j in seq_along(cr.cv)){
	lines(x=cr.mean, y=res[,j], typ='o',col=j)
	text(x=cr.mean[1], y=res[1,j],labels = cr.cv[j],pos=4,col=j)
}

image(x=cr.mean, y=cr.cv, z=res, 
	  col = gray.colors(16),
	  main = 'R', 
	  las=1, xlab='mean',ylab='cv')
contour(x=cr.mean, y=cr.cv, z=res, add = TRUE,nlevels = 6)

best.col <- 'red'
points(x=cr.mean.best, y=cr.cv.best, cex=2, lwd=6, col=best.col)
abline(v=cr.mean.best, col=best.col, lty=4)
abline(h=cr.cv.best, col=best.col, lty=4)
text(x=cr.mean.best, y=cr.cv.best,
	 labels = paste0(round(res[idx],3),' (',round(cr.mean.best,2),';',round(cr.cv.best,2),')'),
	 pos = 3, cex=1.3, col='red')
grid()

dev.off()

t1 <- as.numeric(Sys.time())
msgt <- paste('------ Fit finished in',round((t1-t0)/60,1),'min')
print(msgt)
message(msgt)


