library(plyr)
library(rgexf)
library(igraph)


### ===== Workplace ===== 

# Downloaded from :  http://www.sociopatterns.org/datasets/contacts-in-a-workplace/
x <- read.table('tij_InVS.dat')
x$time <- (x$V1-min(x$V1))/86400
x$hour <- round((x$V1-min(x$V1))/60/60)
x$day  <- round((x$V1-min(x$V1))/60/60/24)
x$key  <- paste(x$V2,x$hour,sep='_')

z  <- subset(x,subset = V2==492)
y  <- ddply(x,c('V2','day', 'key'),summarise, n=length(key))
cr <- ddply(y,c('V2','day'), summarize, cd = length(key))

hist(cr$cd, breaks=0:15,
	 col = 'lightgrey',
	 main = 'Distribution of contact rate at workplace',
	 xlab = 'contacts / person / day')
abline(v=mean(cr$cd), lwd=2, lty=2)
summary(cr$cd)
sd(cr$cd)


###  ==== Primary School =====

summary.contact <- function(file, do.plot = TRUE) {
	
	a <- igraph::read_graph(file = file, format = 'gml')
	dd.x <- igraph::degree(a)
	deg.max <- max(dd.x)
	deg.vec <- 0:deg.max
	dd.y <- igraph::degree_distribution(graph = a)
	mean.deg <- sum( dd.y * deg.vec)
	stddev.deg <- sqrt( sum( dd.y * deg.vec^2) - mean.deg^2 )
	if(do.plot) {
		plot(x=deg.vec, y=dd.y, typ='h',lwd=6, 
			 main = 'Distribution of number of contacts in French primary school',
			 xlab = 'contacts/child/day')
		abline(v=mean.deg, lty=2, lwd=2, col='red')
	}
	return(list(mean.deg=mean.deg, stddev.deg=stddev.deg))
}

# fname <- 'journal.pone.0023176.s003.GML'
par(mfrow=c(1,1))
summary.contact(file= 'journal.pone.0023176.s003.GML')
# SAME???! summary.contact(file= 'journal.pone.0023176.s004.GML')


###  ==== High School =====

# Downloaded from: http://www.sociopatterns.org/datasets/high-school-dynamic-contact-networks/
hs0 <- read.table('thiers_2011.csv')

hs <- data.frame(time = hs0[,1], a=hs0[,2], b=hs0[,3])

hs$time <- hs$time - min(hs$time)

hs$hour <- round(hs$time/60/60)
hs$day  <- round(hs$time/60/60/24)
hs$key  <- paste(hs$a,x$hour,sep='_')

y  <- ddply(hs, c('a','day', 'key'),summarise, n=length(key))
cr.hs <- ddply(y,c('a','day'), summarize, cd = length(key))

hist(cr.hs$cd, breaks=0:max(cr.hs$cd+1),
	 col = 'lightgrey',
	 main = 'Distribution of contact rate at French High school',
	 xlab = 'contacts / person / day')
abline(v=mean(cr.hs$cd), lwd=2, lty=2)
summary(cr.hs$cd)
sd(cr.hs$cd)
