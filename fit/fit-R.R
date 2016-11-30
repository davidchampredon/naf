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
n.MC  <- 2
n.cpu <- 2

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
cr.mean <- seq(1.0, 3.0, by=1)
cr.sd   <- seq(0.5, 3.0, by=1)

n.cpu.cr.mean <- 2

res <- matrix(nrow = length(cr.mean), ncol=length(cr.sd))
n <- prod(dim(res))

# ==== Simulation Loop ====
# Run the simulations on all scenarios:

wrap_sim <- function (i, cr.mean, cr.sd, 
					  simul.prm, prm, n.MC, n.cpu)
{
	for(j in seq_along(cr.sd)){
		# Overwrite parameter values
		# associated with current scenario:
		prm[['contact_rate_mean']]   <- cr.mean[i]
		prm[['contact_rate_stddev']] <- cr.sd[j]
		
		# Run the simulation for that scenario:
		stoch_build_world <- simul.prm[['build_world_stoch']]
		res.list.0 <- run.simul.fit.R(n.MC,n.cpu)	
		
		# merge all MC iterations:
		pop   <- merge.pop.mc(res.list = res.list.0,
							  n.cpu = n.cpu, 
							  doparallel = TRUE)
		
		# Filter out fizzles:
		idx.mc         <- unique(pop$mc)
		fizz.mc        <- identify.fizzle(pop)
		fizz.mc.idx    <- which(fizz.mc==TRUE)
		if(length(fizz.mc.idx)==0) idx.mc.no.fizz <- idx.mc
		if(length(fizz.mc.idx)>0) idx.mc.no.fizz <- idx.mc[-fizz.mc.idx]
		
		if(length(idx.mc.no.fizz) > 0){
			R <- calc.R(pop)
			res[i,j] <- R[['R.mean']]
		}
		if(length(idx.mc.no.fizz) == 0) res[i,j] <- NA
		return(res[i,])
	}
}

# a <- wrap_sim(i=3, cr.mean, cr.sd, simul.prm, prm, n.MC, n.cpu)

# Parallel computation:
sfInit(parallel = TRUE , cpu = n.cpu.cr.mean)
sfLibrary(naf, lib.loc = R.library.dir)
sfLibrary(snowfall)
sfLibrary(plyr)
sfExportAll()
res.tmp <- sfSapply(x = seq_along(cr.mean), fun = wrap_sim, 
					cr.mean=cr.mean, cr.sd=cr.sd,
					simul.prm=simul.prm,
					prm=prm, n.MC=n.MC, 
					n.cpu=n.cpu,
					simplify = FALSE)
sfStop()

res <- matrix(unlist(res.tmp),nrow = length(cr.mean))

save.image(file = 'tmp.RData')

# ==== Best Fit ====

dist <- abs(res - R.target)
idx <- which(dist == min(dist,na.rm = TRUE), arr.ind = TRUE)
res[idx]
cr.mean.best <- cr.mean[idx[1]]
cr.sd.best   <- cr.sd[idx[2]]

write.csv(x = data.frame(cr.mean.best=cr.mean.best, cr.sd.best=cr.sd.best),
		  file = 'fitted_cr.csv', quote = F, row.names = F)

# ==== PLOTS ====

pdf('fitted_cr.pdf')

plot(1,xlim=range(cr.sd), ylim=range(res,na.rm = TRUE),
	 xlab = 'sd', ylab='R')
for(i in seq_along(cr.mean)){
	lines(x=cr.sd, y=res[i,], typ='o', col=i)
	text(x=cr.sd[1], y=res[i,1],labels = cr.mean[i],pos=4,col=i)
}

plot(1,xlim=range(cr.mean), ylim=range(res,na.rm = TRUE),
	 xlab = 'mean', ylab='R')
for(j in seq_along(cr.sd)){
	lines(x=cr.mean, y=res[,j], typ='o',col=j)
	text(x=cr.mean[1], y=res[1,j],labels = cr.sd[j],pos=4,col=j)
}

image(x=cr.mean, y=cr.sd, z=res, 
	  col = heat.colors(12),
	  main = 'R', 
	  las=1, xlab='mean',ylab='sd')
contour(x=cr.mean, y=cr.sd, z=res, add = TRUE,nlevels = 6)
points(x=cr.mean.best, y=cr.sd.best, cex=3, lwd=3)
abline(v=cr.mean.best)
abline(h=cr.sd.best)
grid()
dev.off()

t1 <- as.numeric(Sys.time())
msgt <- paste('------ Fit finished in',round((t1-t0)/60,1),'min')
print(msgt)


