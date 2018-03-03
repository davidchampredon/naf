#############################################
#####
#####  FRAILTY SHAPE PARAMETERS CALIBRATION (MANUALLY!)
#####
#############################################

# Retrieve Data 
dat <- read.csv('../data/frailty/frailty-canada.csv')

# Shape parameters
F0   <- 0.60   # Frailty at age 0
q    <- 1.0    # Shape first part of the frailty curve # 2.0
amin <- 30     # Age where frailty is minimal #30
Fmin <- 0.12   # Minimum frailty value
ap   <- 60     # age at inflexion point for second part of the frailty curve
b    <- 1.3    # Shape second part of the frailty curve

age <- seq(0,100,by=0.1)

frailty.age <- function(age,F0,q,amin,Fmin,ap,b) {
	
	age.infant <- 2
	yy0 <- (F0-Fmin)/age.infant^q * (age.infant-age)^q + Fmin
	yy0[age>age.infant] <- 0
	# plot(yy0)
	
	yy1 <- yy0
	yy1[age <= age.infant] <- 0
	yy1[age>age.infant & age<amin] <- Fmin
	# plot(yy1)
	
	alpha <- -log(1/Fmin/b-1/b)/(amin-ap)
	yy2 <- 1/(1+b*exp(-alpha*(age-ap)))
	yy2[age<amin] <- 0
	f <- yy0 + yy1 + yy2
	return(f)
}

fr0 <- frailty.age(age, F0, q,
				   amin,
				   Fmin,
				   ap,
				   b) 
# par(mfrow=c(1,1)); plot(age,fr0,typ='l')

	
	
dist.dat <- function(x,F0,q,dat) {
	fr <- frailty.age(age=dat$age, F0=F0,q=q,
					  amin=x[['amin']],Fmin=x[['Fmin']],ap=x[['ap']],b=x[['b']])
	dist <- 0
	for(i in 1:nrow(dat)) dist <- dist + (fr[i]-dat$frailty)^2
	return(sum(dist))
}
x0 <- c(amin=amin,Fmin=Fmin,ap=ap,b=b)
dist.dat(x0, F0,q,dat)

fit <- optim(par = x0, 
			 fn = dist.dat, #method = 'SANN',
			 F0=F0, q=q, dat=dat)

fit.prm <- fit$par
print(fit.prm)

# ==== Plot ====

fr0 <- frailty.age(age, F0, q,
				  amin,
				  Fmin,
				  ap,
				  b) 

fr <- frailty.age(age, F0, q,
				  amin = fit.prm[['amin']],
				  Fmin = fit.prm[['Fmin']],
				  ap   = fit.prm[['ap']],
				  b    = fit.prm[['b']]) 

myplot <- function(age, fr, dat, title) {
	plot(age,fr,typ='l', ylim=c(0,1), 
		 xlab = 'age', ylab = 'frailty',
		 main = title,
		 las=1,lwd=6)
	# data:
	points(dat$age, dat$frailty, pch=1, lwd=2, col='red')
	
	abline(v=amin,lty=2,col='grey',lwd=1)
	abline(h=Fmin,lty=2,col='grey',lwd=1)
	abline(v=ap,lty=4,col='orange',lwd=1)
	yap <- fr[age==ap]
	abline(h=yap,lty=4,col='orange',lwd=1)
	grid()
}

par(mfrow=c(1,2))
myplot(age, fr, dat, title='Fit with optim')
myplot(age, fr0, dat, title='Manual fit')



