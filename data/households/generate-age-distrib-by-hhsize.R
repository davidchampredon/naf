###
###   GENERATE SYNTHETIC AGE DISTRIBUTIONS 
###   CONDITIONAL ON HOUSEHOLD SIZE
###

remove.Inf <- function(x, mult = 2) {
	z <- x
	zz <- z
	z[z==Inf] <- 0
	zz[zz==Inf]	 <- mult * max(z)
	return(zz)
}


gen.ad <- function(prm, do.plot = TRUE, plot.title=''){
	
	agemean <- prm[['agemean']]
	agemin  <- prm[['agemin']]
	agemax  <- prm[['agemax']]
	a       <- prm[['a']]
	
	qv <- vector()
	pv <- c(0.05, 0.50,0.95,0.99,0.999)
	
	xbar <- (agemean-agemin)/(agemax-agemin)
	b <- a/xbar -a
	x = seq(0,1,by=1/(agemax-agemin))
	age <- agemin + (agemax-agemin)*x
	ad <- dbeta(x, shape1 = a, shape2 = b)
	ad <- remove.Inf(ad)
	
	if(do.plot){
		plot(age, ad, typ='l',lwd=6, xlim=c(0,100),
			 main = plot.title,
			 cex.axis=0.7,
			 cex.main = 0.8,
			 mgp=c(1,0.1,0),
			 yaxt='n',ylab='',xlab='age')
		
		qv <- agemin + (agemax-agemin) * qbeta(p=pv, shape1 = a, shape2 = b)
		abline(v=qv, lty=2, col='tomato', lwd=2)
		text(x=qv,y=0.91^(1:length(qv))*max(ad),
			 cex = 0.7,
			 labels = paste0(round(qv),'(',pv*100,'%)'),
			 pos = 4, col='tomato')
		#grid()
	}
	return(list(age=age, ad=ad))
}



# FOR QUICK TEST:
if(0){ 
	# agemin  agemax agemean       a 
	# 18.0    65.0    66.0     0.7 
	
	par(mfrow=c(1,1))
	aa <- gen.ad(prm=c(agemin=18, agemax=65, agemean=66,   a=0.7),
				 do.plot = TRUE, plot.title='test')

}



M <- list()
M[[1]] <- list(c(agemin=18,agemax=80, agemean=60,  a=1.5))

M[[2]] <- list(c(agemin=18,agemax=80, agemean=50,  a=2.5),
			   c(agemin=1 ,agemax=80, agemean=50,  a=2.5))

M[[3]] <- list(c(agemin=18,agemax=80, agemean=35,  a=1.8),
			   c(agemin=18,agemax=80, agemean=35,  a=1.8),
			   c(agemin=1, agemax=50, agemean=10,  a=2.0))

shift.mean <-  c(agemin=0, agemax=0, agemean= 2.9,   a=0)

M[[4]] <- list(M[[3]][[1]] + shift.mean,
			   M[[3]][[2]] + shift.mean,
			   M[[3]][[3]] + shift.mean,
			   M[[3]][[3]])

M[[5]] <- list(M[[4]][[1]] + shift.mean,
			   M[[4]][[2]] + shift.mean,
			   M[[4]][[3]] + shift.mean,
			   M[[4]][[4]] + shift.mean,
			   M[[4]][[4]])

M[[6]] <- list(M[[5]][[1]] + shift.mean,
			   M[[5]][[2]] + shift.mean,
			   M[[5]][[3]] + shift.mean,
			   M[[5]][[4]] + shift.mean,
			   M[[5]][[5]] + shift.mean,
			   M[[5]][[5]])



save.ad <- function(M, hh, indiv, fname,...) {
	z <- gen.ad(prm = M[[hh]][[indiv]] )#, ...)
	zz <- matrix(unlist(z),ncol = 2)
	zz <- zz[!is.infinite(zz[,2]),]
	zz[,2] <- zz[,2]/sum(zz[,2])
	write.table(zz,
				file=paste0(paste(fname,hh-1,indiv-1,sep="_"),'.csv'), # '-1' to be consistent with C++ code
				sep = ',',
				col.names = FALSE,
				row.names = FALSE)
}

save.ad.to.file <- TRUE

if(save.ad.to.file) pdf('plot_age_distrib_hhsize.pdf', width=20,height = 20)
par(mfrow=c(6,6), mar=rep(0.9,4))
for(i in 1:6){
	for(j in 1:6){
		if (j < i) plot.new()
		if(j >= i) {
			ptt <- paste0('HH size=',j,' ; indiv #',i)
			#gen.ad(M[[j]][[i]],plot.title = ptt)
			save.ad(M, j, i, 'hh_size_ad', plot.title=ptt)
		}
	}
}
if (save.ad.to.file) dev.off()






