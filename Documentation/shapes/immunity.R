
# Humoral immunity

q <- 100
x <- seq(0,q,length.out = 1e3)

plot.hum <- function(x,q,p,h0, add) {
	y <- h0 * ( (q^p-(x)^p)^(1/p) )/q
	if(!add){
		plot(x,y,typ='l',lwd=6, 
			 main='Humoral immunity',
			 ylab='',xlab='age',
			 ylim = c(0,1),las=1)
		grid(lty=2)	}
	if(add){
		lines(x,y,typ='l',lwd=2)
	}
}

# Cellular immunity
plot.cel <- function(imm.max, slope, pivot,add) {
	agemax <- 100
	x <- seq(0,agemax,length.out = 1000)
	y <- imm.max/(1+exp(-slope*(x/pivot-1)))
	
	if(!add){
		plot(x,y,typ='l', xlim=c(0,agemax), 
			 main='Cellular immunity index',
			 las=1, 
			 xlab='Age (years)',ylab='',
			 ylim=c(0,1),lwd=6)
		grid(lty=2)
	}
	if(add){
		lines(x,y,lwd=2)
	}
}

save.plot <- TRUE
if(save.plot) pdf('../figures/immunity.pdf',width = 12, height = 6)
par(mfrow=c(1,1), cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
# h0 <- 0.5
# plot.hum(x,q=100,p=1.7,h0, add=F)
# plot.hum(x,q=100,p=3.0,  h0, add=T)
plot.cel(imm.max = 0.39, slope=2, pivot=20, add=F)
# plot.cel(imm.max = 0.7, slope=2, pivot=29, add=T)
if(save.plot) dev.off()


