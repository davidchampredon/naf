x <- y <- seq(0,100,by=1)


f <- function(x,y,p,q,r) {
	exp(-p*(x-y-q)^2)*exp(-r*(x^2+y^2))    #/(1+r*(x^2+y^2))
}

z <- matrix(nrow = length(y), ncol = length(x))

z <- z / sum(z)

tmp <- vector()

for(j in seq_along(x)){
	for(i in seq_along((y))){
		r <- 4*1e-4
		tmp[1] <- f(x[j],y[i],p=1e-2,q=0, r)
		tmp[2] <- f(x[j],y[i],p=3e-2,q=30, r)
		tmp[3] <- f(x[j],y[i],p=3e-2,q=-30, r)
		tmp[4] <- f(x[j],y[i],p=1e-2,q=60, r)
		tmp[5] <- f(x[j],y[i],p=1e-2,q=-60, r)
		z[i,j] <- sum(tmp)
	}
	
}

image(x=x, y=y, z=z, col = topo.colors(12,0.89))
grid(lty=1, col=rgb(1,1,1,0.2))
