###
###  Shape for the probability to die
###

frailty <- seq(0,1,by=0.01)

# Shape parameters:
pmin <- 0.15
pmax <- 0.5
f.cvx <- 0.4
b <- 10

# Death probability:

mult <- pmin + 1/(1+exp(-b*(frailty-f.cvx))) * (pmax-pmin)

p <- mult * frailty

# Plot:
par(mfrow=c(1,2))
plot(x=frailty, y=mult, 
	 main = 'Death probability multiplier',
	 ylab='',
	 typ = 'l',	 lwd = 6, ylim=c(0,1), 
	 las = 1)
grid()
abline(h=c(pmin,pmax), lty =2)
abline(v=f.cvx, lty =2)

plot(x=frailty, y=p, 
	 main = 'Death probability',
	 ylab='',
	 typ = 'l',	 lwd = 6, ylim=c(0,1), 
	 las = 1)
abline(a=0,b=1, lty=2)
abline(a=0,b=pmax, lty=3)
grid()