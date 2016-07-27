z <- 0.1
a <- 1/z
b = 2/z
x <- seq(0,1,length.out = 1000)

y <- 1/(1+exp(-(b*x-a)))


plot(x,y,typ='l', xlim=c(0,1), ylim=c(0,1),lwd=3)
grid(lty=2)