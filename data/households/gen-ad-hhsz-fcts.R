###
###   GENERATE SYNTHETIC AGE DISTRIBUTIONS 
###   CONDITIONAL ON HOUSEHOLD SIZE
###


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

### Generate all age distributions given household size.
### 'agemean.vec' and 'a.vec' define the beta dstributions.
gen.all.ad.hhsz <- function(agemean.vec, a.vec, path.save='./') {
	M <- list()
	M[[1]] <- list(c(agemin=18,agemax=80, agemean=agemean.vec[1],  a=a.vec[1]))
	
	M[[2]] <- list(c(agemin=18,agemax=80, agemean=agemean.vec[2],  a=a.vec[2]),
				   c(agemin=18,agemax=80, agemean=agemean.vec[3],  a=a.vec[3]))
	
	M[[3]] <- list(c(agemin=18,agemax=65, agemean=agemean.vec[4],  a=a.vec[4]),
				   c(agemin=18,agemax=65, agemean=agemean.vec[5],  a=a.vec[5]),
				   c(agemin=0, agemax=40, agemean=agemean.vec[6],  a=a.vec[6]))
	
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
		z <- gen.ad(prm = M[[hh]][[indiv]], ...)
		zz <- matrix(unlist(z),ncol = 2)
		zz <- zz[!is.infinite(zz[,2]),]
		zz[,2] <- zz[,2]/sum(zz[,2])
		write.table(zz,
					file=paste0(paste(fname,hh-1,indiv-1,sep="_"),'.csv'), # '-1' to be consistent with C++ code
					sep = ',',
					col.names = FALSE,
					row.names = FALSE)
	}
	for(i in 1:6){
		for(j in 1:6){
			if(j >= i) {
				ptt <- paste0('HH size=',j,' ; indiv #',i)
				save.ad(M, j, i, 
						fname = paste0(path.save,'hh_size_ad'), 
						do.plot = FALSE,
						plot.title=ptt)
			}
		}
	}
	
}



if(F){ ### FOR TESTING:
	
	agemean.vec <- c(60,50,50,35,35,10)
	a.vec       <- c(1.5,2.5,2.5,1.8,1.8,2)
	
	par(mfrow=c(4,6), mar=rep(0.9,4))
	gen.all.ad.hhsz(agemean.vec, a.vec)
}


read.all.ad.hhsz <- function(path) {
	x <- system(paste0('ls ',path,'hh_size_ad*.csv'), intern = TRUE)
	
	par(mfrow=c(4,6), mar=rep(0.9,4))
	for (i in seq_along(x)){
		a <- read.csv(file = x[i], header = FALSE)
		
		qv <- vector()
		pv <- c(0.05, 0.50,0.95,0.99,0.999)
		title <- gsub(path,'',x[i])
			plot(a[,1], a[,2], typ='l',lwd=6, xlim=c(0,100),
				 main = title,
				 cex.axis=0.7,
				 cex.main = 0.8,
				 mgp=c(1,0.1,0),
				 yaxt='n',ylab='',xlab='age')
			grid()
	}
}




