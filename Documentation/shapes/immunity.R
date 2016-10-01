
# Humoral immunity

q <- 100
p <- 2
x <- seq(0,q,length.out = 1e3)


plot.hum <- function(x,q,p, add) {
	y <- 1 * ( (q^p-(x)^p)^(1/p) )/q
	if(!add){
		plot(x,y,typ='l',lwd=6, 
			 main='Humoral immunity',
			 ylab='',xlab='age',
			 ylim = c(0,1),las=1)
		grid(lty=1)	}
	if(add){
		lines(x,y,typ='l',lwd=2)
	}
}

# Cellular immunity
plot.cel <- function(slope,add) {
	agemax <- 100
	# slope <- 3	
	pivot <- 20
	b <- slope/pivot
	x <- seq(0,agemax,length.out = 1000)
	imm.max <- 0.7
	y <- imm.max/(1+exp(-slope*(x/pivot-1)))
	
	if(!add){
		plot(x,y,typ='l', xlim=c(0,agemax), 
			 main='Cellular immunity',
			 las=1, 
			 xlab='age',ylab='',
			 ylim=c(0,1),lwd=6)
		grid(lty=2)
	}
	if(add){
		lines(x,y,lwd=2)
	}
}


pdf('../figures/immunity.pdf',width = 12, height = 6)
par(mfrow=c(1,2))
plot.hum(x,q=100,p=2, add=F)
plot.hum(x,q=100,p=3, add=T)
plot.cel(slope=3, add=F)
plot.cel(slope=1.5, add=T)
dev.off()


