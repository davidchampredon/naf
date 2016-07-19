# Simple simulation:
prm1 <- c(
	infectious_mean = 3,
	popSize         = 4000,
	R0              = 3)

prm2 <-  c(
	horizon = 150,
	latent_mean     = 2,	
	nE = 1,
	nI = 1,
	init_I1 = 1,
	n.time.steps = 500,
	per.capita = FALSE
)

# Synthetic data
sim <- simul.SEmInR.det(prm1, prm2)
df  <- sim$ts
det.inc <- df$inc
det.prev <- df$Iall
det.t <- df$time

plot(det.t,det.prev,typ='l')
lines(ts$time, ts$prevalence, col='red')
