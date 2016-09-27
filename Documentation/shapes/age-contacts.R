x <- y <- seq(0,100,by=0.5)


psi <- function(x,y,a,b,w,q,r) {
	exp( -( ( (x-a) -(y-b) - q )^2  )/w^2 - ((x-a)^2+(y-b)^2)/r^2) 
}

z <- matrix(nrow = length(y), ncol = length(x))
tmp <- vector()

for(j in seq_along(x)){
	for(i in seq_along((y))){
		tmp[1] <- psi(x[j],y[i],a=0,b=0,w=9,q=28, r=45)   # parents/children
		tmp[2] <- psi(x[j],y[i],a=0,b=0,w=9,q=-28,r=45)   # parents/children
		
		tmp[3] <- psi(x[j],y[i],a=10,b=10,w=9,q=0,r=28)    # children/children
		
		tmp[4] <- psi(x[j],y[i],a=45,b=45,w=20,q=0,r=17)    # adults/adults
		
		tmp[5] <- psi(x[j],y[i],a=75,b=75,w=15,q=0,r=15)    # seniors
		
		tmp[6] <- psi(x[j],y[i],a=0,b=0,w=10,q=60,r=55)    # seniors/children
		tmp[7] <- psi(x[j],y[i],a=0,b=0,w=10,q=-60,r=55)    # seniors/children
		
		z[i,j] <- sum(tmp)
	}
}
z <- z / sum(z)

# PLOTS

par(mfrow=c(2,2))

image(x=x, y=y, z=z, 
	  col = topo.colors(12,0.5),
	  main = 'Age contact preference matrix', 
	  las=1, xlab='age',ylab='')
contour(x=x, y=y, z=z, add = TRUE,nlevels = 6)
grid()
abline(a=0,b=1,lty=2)

col.b <- rgb(0,0,0,0.1)
col.f <- rgb(0,0,0,0.1)

persp(x=x, y=y, z=z,
	  phi = 10, 
	  theta = -45,
	  col = col.f, border = col.b)

persp(x=x, y=y, z=z,
	  phi = 50, 
	  theta = 0,
	  col = col.f, border = col.b)
persp(x=x, y=y, z=z,
	  phi = 10, 
	  theta = 45,
	  col = col.f, border = col.b)
