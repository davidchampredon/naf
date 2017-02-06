
# ---- Humoral immunity ----  

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

# ---- Cellular immunity ----

# Data retrieved from Fig 4 
# Beest te DE, Birrell PJ, Wallinga J, De Angelis D, van Boven M. Joint modelling of serological and hospitalization data reveals that high levels of pre-existing immunity and school holidays shaped the influenza A pandemic of 2009 in The Netherlands. Journal of The Royal Society Interface 2014; 12: 20141244–4.
M <- rbind(c(0,4,0,0),
           c(5,9, 0, 0.3),
           c(10,19, 0.2, 0.45),
           c(20,64, 0.6, 0.8),
           c(65,100, 0.6, 0.95))
colnames(M) <- c('x0','x1','y0','y1')
M <- as.data.frame(M)
#Median value:
M.med.x <- c(2.5, 7, 15, 42, 80)
M.med.y <- c(0, 0.09, 0.29, 0.72, 0.92)

imm.cell <- function(age, slope, imm.max, pivot) {
    return(imm.max/(1+exp(-slope*(age/pivot-1))))
}

init <- c(slope=2, pivot=25, imm.max=0.9)

ferr <- function(x){
    valf <- imm.cell(age=M.med.x, slope=x[['slope']], 
                     imm.max = x[['imm.max']], pivot=x[['pivot']])
    return(sum((valf-M.med.y)^2))
}

xstar <- optim(par = init, fn = ferr)$par


plot.cel <- function(imm.max, slope, pivot,add) {
	agemax <- 100
	x <- seq(0,agemax,length.out = 1000)
	y <- imm.cell(age=x, slope, imm.max, pivot)
	
	if(!add){
		
	    plot(x,y,typ='l', xlim=c(0,agemax), 
			 main='',
			 las=1, 
			 xlab='Age (years)',ylab='Cellular immunity level',
			 ylim=c(0,1),lwd=6)
	    
	   for(i in 1:nrow(M)){
	       xvec <- c(M$x0[i], M$x1[i])
	       polygon(x = c(xvec,rev(xvec)), 
	               y = c(M$y0[i],M$y0[i],M$y1[i],M$y1[i]),
	               col = rgb(0,0,0,0.05), border = NA)
	   }
	    
	   points(x=M.med.x, y=M.med.y, pch=1, cex=2, lwd=3)
		#grid(lty=2)
	}
	if(add){
		lines(x,y,lwd=2)
	}
}

save.plot <- T
if(save.plot) pdf('../figures/immunity.pdf',width = 12, height = 6)
par(mfrow=c(1,1), cex.lab=1.1, cex.axis=1.5, cex.main=1.5)
# h0 <- 0.5
# plot.hum(x,q=100,p=1.7,h0, add=F)
# plot.hum(x,q=100,p=3.0,  h0, add=T)
plot.cel(imm.max = xstar[['imm.max']], 
         slope   = xstar[['slope']], 
         pivot   = xstar[['pivot']], 
         add=F)
xstar
# plot.cel(imm.max = 0.7, slope=2, pivot=29, add=T)
if(save.plot) dev.off()


