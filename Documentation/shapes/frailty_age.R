
F0 <- 0.6
ap<- 60
amin <- 30
Fmin <- 0.15
b <- 2
alpha <- -log(1/Fmin/b-1/b)/(amin-ap)

xx <- seq(0,100,by=0.1)
yy <- 1/(1+b*exp(-alpha*(xx-ap)))
yy[xx<=amin] <- NA

q <- 3
yy0 <- (F0-Fmin)/amin^q * (amin-xx)^q + Fmin
yy0[xx>amin] <- NA

plot(xx,yy,typ='l', ylim=c(0,1), las=1,lwd=6)
lines(xx,yy0, lwd=6)

abline(v=amin,lty=2,col='grey',lwd=2)
abline(h=Fmin,lty=2,col='grey',lwd=2)
abline(v=ap,lty=4,col='orange',lwd=2)
yap <- yy[xx==ap]
abline(h=yap,lty=4,col='orange',lwd=2)
grid()


