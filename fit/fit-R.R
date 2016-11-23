###
###   FIT CONTACT RATE TO R0
###

t0 <- as.numeric(Sys.time())

# Load all default parameters:
source('fit_fcts.R')
source('../simul/utils-run.R')
source('../simul/utils-misc.R')

R.library.dir   <- '../Rlibrary/lib'
param.model.dir <- '../param-model/'
data.dir        <- '../data/'
simul.dir       <- '../simul/'

PRM <- load.all.parameters(R.library.dir ,
						   param.model.dir,
						   data.dir,
						   simul.dir) 
for(i in seq_along(PRM)){
	assign(x = names(PRM)[i], value = PRM[[i]])
}

n.MC  <- 10
n.cpu <- 10

cr.mean <- seq(1.0, 4.0, by=1)
cr.sd   <- seq(0.5, 4.0, by=0.5)

res <- matrix(nrow = length(cr.mean), ncol=length(cr.sd))

# Run the simulations on all scenarios:
for(i in seq_along(cr.mean)){
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
		R <- calc.R(pop)
		res[i,j] <- R[['R.mean']]
	}
}
save.image(file = 'tmp.RData')

plot(1,xlim=range(cr.mean), ylim=range(res))
for(j in seq_along(cr.sd)){
	lines(x=cr.mean, y=res[,j], typ='o')
	text(x=cr.mean[1], y=res[1,j],labels = cr.sd[j],pos=2)
}


image(x=cr.mean, y=cr.sd, z=res, 
	  col = heat.colors(12),
	  main = 'R', 
	  las=1, xlab='mean',ylab='sd')
contour(x=cr.mean, y=cr.sd, z=res, add = TRUE,nlevels = 6)
grid()

t1 <- as.numeric(Sys.time())
msgt <- paste('------ Fit finished in',round((t1-t0)/60,1),'min')
print(msgt)






