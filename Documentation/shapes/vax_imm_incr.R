save.to.file <- T
m <- 0.70
v <- m/100

x<-seq(0,2,length.out = 1000)

y <- dgamma(x,shape = m^2/v, rate = m/v)
q <- qgamma(p = c(0.025, 0.975),shape = m^2/v, rate = m/v)

if(save.to.file) pdf('../figures/vax_imm_incr.pdf')
plot(pmin(1,x),y,typ='l', lwd=6, 
	 xlab = 'Immunity increase', ylab='', yaxt='n',
	 main ='Distribution shape of immunity \n increase following vaccination')
abline(v=c(m,q),lty=2)
grid()
if(save.to.file) dev.off()
