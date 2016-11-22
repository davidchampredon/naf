###
###  Contact rate before transmission
###

set.seed(1234)

# Contact rate

cr_mean <- 7
cr_sd   <- 5.0

tmp  <- 1 + cr_sd*cr_sd/cr_mean/cr_mean
m_ln <- log(cr_mean/sqrt(tmp))
s_ln <- sqrt(log(tmp))

x <- seq(0,30,length.out = 1e3)
y <- dlnorm(x,meanlog = m_ln,sdlog = s_ln)

pq <- c(0.01,0.25,0.5,0.75,0.99)
q <- qlnorm(p=pq,meanlog = m_ln,sdlog = s_ln)

# Number of contacts:
n <- 1e5
cr.sample <- rlnorm(n=n, meanlog = m_ln,sdlog = s_ln)
num.c <- rpois(n=n, lambda = cr.sample)
qpoiss <- quantile(num.c,probs = pq)
	
# Plots

#pdf('contact-rate.pdf',width = 10)
par(mfrow=c(1,2))
plot(x,y,typ='l',
	 main = 'Contact Rate Distribution',
	 lwd=6,yaxt='n', bty="n",
	 xlab='contacts/day/capita',ylab='')
abline(v=q, col='gold',lty=2)

hist(num.c,breaks = 0:max(num.c),
	 yaxt='n', ylab='', xlab='contacts/day/capita',
	 main = 'Number of contacts',
	 border = 'lightgrey', col='gray')
abline(v=qpoiss, col='gold',lty=2)
#dev.off()


