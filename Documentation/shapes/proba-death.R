###
###  Shape for the probability to die
###

frailty <- seq(0,1,by=0.01)

# Shape parameters:
pmin <- 0.15
pmax <- 0.95
f.cvx <- 0.4
b <- 10

# Death probability:
p <- pmin + 1/(1+exp(-b*(frailty-f.cvx))) * (pmax-pmin)

# Plot:
plot(x=frailty, y=p, 
	 main = 'Death probability multiplier',
	 ylab='',
	 typ = 'l',	 lwd = 6, ylim=c(0,1), 
	 las = 1)
grid()