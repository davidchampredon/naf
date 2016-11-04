#############################################
#####
#####  FRAILTY SHAPE PARAMETERS CALIBRATION (MANUALLY!)
#####
#############################################



# Retrieve Data 
dat <- read.csv('../data/frailty/frailty-canada.csv')

# Shape parameters
F0 <- 0.25
ap<- 60
amin <- 30
Fmin <- 0.12
b <- 1.3
alpha <- -log(1/Fmin/b-1/b)/(amin-ap)

xx <- seq(0,100,by=0.1)
yy <- 1/(1+b*exp(-alpha*(xx-ap)))
yy[xx<=amin] <- NA

q <- 3
yy0 <- (F0-Fmin)/amin^q * (amin-xx)^q + Fmin
yy0[xx>amin] <- NA

# ==== Plot ====

plot(xx,yy,typ='l', ylim=c(0,1), 
	 xlab = 'age', ylab = 'frailty',
	 las=1,lwd=6)
lines(xx,yy0, lwd=6)
# data:
points(dat$age, dat$frailty, pch=1, lwd=2, col='red')

abline(v=amin,lty=2,col='grey',lwd=1)
abline(h=Fmin,lty=2,col='grey',lwd=1)
abline(v=ap,lty=4,col='orange',lwd=1)
yap <- yy[xx==ap]
abline(h=yap,lty=4,col='orange',lwd=1)
grid()

