###
###  Shape for the probability to die
###

frailty <- seq(0,1,by=0.01)

# Shape parameters:
pmin <- 0.1
pmax <- 0.9
f.cvx <- 0.6
b <- 15

# Death probability:
p <- pmin + 1/(1+exp(-b*(frailty-f.cvx))) * (pmax-pmin)

# Plot:
plot(x=frailty, y=p, 
	 main = 'Death probability multiplier',
	 ylab='',
	 typ = 'l',	 lwd = 6, ylim=c(0,1), 
	 las = 1)
grid()