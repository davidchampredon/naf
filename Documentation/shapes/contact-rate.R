###
###  Contact rate before transmission
###

set.seed(1234)

# Contact rate

mult    <- 1
cr_mean <- 6 * mult
cr_sd   <- 3

# Lognormal parameterization
tmp  <- 1 + cr_sd*cr_sd/cr_mean/cr_mean
m_ln <- log(cr_mean/sqrt(tmp))
s_ln <- sqrt(log(tmp))

# Gamma parameterization
shape <- cr_mean^2 / cr_sd^2
scale <- cr_sd^2 / cr_mean

x <- seq(0,30,length.out = 1e3)
y <- dlnorm(x,meanlog = m_ln,sdlog = s_ln)
y.g <- dgamma(x,shape = shape, scale = scale)

pq <- c(0.01,0.5,0.99)
q   <- qlnorm(p = pq, meanlog = m_ln,sdlog = s_ln)
q.g <- qgamma(p = pq, shape = shape, scale = scale)

# Number of contacts:
n <- 1e5
cr.sample.ln <- rlnorm(n=n, meanlog = m_ln,sdlog = s_ln)
cr.sample.ga <- rgamma(n=n, shape = shape, scale = scale)

num.c.ln <- rpois(n=n, lambda = cr.sample.ln)
num.c.ga <- rpois(n=n, lambda = cr.sample.ga)

qpoiss.ln <- quantile(num.c.ln, probs = pq)
qpoiss.ga <- quantile(num.c.ga, probs = pq)
	
# Plots

#pdf('contact-rate.pdf',width = 10)
# par(mfrow=c(1,2))
layout(matrix(c(1,3,2,3), 2, 2, byrow = TRUE))

plot.density <- function(log) {

	plot(x,y,typ='l',
		 main = 'Contact Rate Distribution',
		 lwd=3,yaxt='n', bty="n",
		 xlab='contacts/day/capita',ylab='',
		 log=log)
	abline(v=q, col='black',lty=2)
	
	lines(x,y.g, col='red', lwd=3)
	abline(v=q.g, col='red',lty=2)
	
	legend(x='topright',legend = c('lognorm','gamma'),col = c('black','red'), lwd=2)
}

plot.density(log='')
plot.density(log='y')

hist(num.c.ln,breaks = 0:max(num.c.ln), 
	 yaxt='n', ylab='', xlab='contacts/day/capita',
	 main = 'Number of contacts',
	 border = 'lightgrey', col=rgb(0,0,0,0.2))
abline(v=qpoiss.ln, col='grey',lty=2)
points(x=max(num.c.ln),y=0,pch=16)

hist(num.c.ga,breaks = 0:max(num.c.ga), add = TRUE,
	 yaxt='n', ylab='', xlab='contacts/day/capita',
	 main = 'Number of contacts',
	 border = 'pink', col=rgb(1,0,0,0.2))
abline(v=qpoiss.ga, col='pink',lty=2)
points(x=max(num.c.ga),y=0,pch=16, col='red')
#dev.off()


