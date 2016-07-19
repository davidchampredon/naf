
library(plyr)
library(gridExtra)

analysis.population <- function(pop) {
	
	sum(pop$is_infected)
	sum(pop$is_infectious)
	sum(pop$is_recovered)
	sum(pop$is_latent)
}

plot.epi.timeseries <- function(ts){
	g <- ggplot(ts, aes(x=time))
	g <- g + geom_line(aes(y=nS),colour='springgreen3') 
	g <- g + geom_line(aes(y=nR),colour='blue')
	
	g1 <- ggplot(ts, aes(x=time))
	g1 <- g1 + geom_line(aes(y=nE),colour='orange') 
	g1 <- g1 + geom_line(aes(y=nIs),colour='red')
	
	g2 <- ggplot(ts, aes(x=time))
	g2 <- g2 + geom_step(aes(y=incidence))
	g2 <- g2 + geom_line(aes(y=prevalence))
	
	grid.arrange(g,g1,g2,nrow=2)
}
# 
# inc <- ts$incidence
# prev <- ts$prevalence
# inc2 <- diff(prev)
# 
# plot(inc, typ='o', xlim=c(0,100))
# lines(c(0,inc2),col='red')
